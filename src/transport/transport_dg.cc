/*!
*
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
* 
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License version 3 as published by the
* Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
* 
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* 
* @file    transport_dg.cc
* @brief   Discontinuous Galerkin method for equation of transport with dispersion.
* @author  Jan Stebel
*/

#include "system/index_types.hh"
#include "system/sys_profiler.hh"
#include "transport/transport_dg.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/dh_cell_accessor.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"
#include "la/linsys_PETSC.hh"
#include "coupling/balance.hh"
#include "coupling/generic_assembly.hh"
#include "transport/advection_diffusion_model.hh"
#include "transport/concentration_model.hh"
#include "transport/heat_model.hh"
#include "transport/assembly_dg.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "input/factory.hh"
#include "fields/equation_output.hh"
#include "mesh/accessors.hh"

FLOW123D_FORCE_LINK_IN_CHILD(concentrationTransportModel)
FLOW123D_FORCE_LINK_IN_CHILD(heatModel)



using namespace Input::Type;

template<class Model>
const Selection & TransportDG<Model>::get_dg_variant_selection_input_type() {
    return Selection("DG_variant", "Type of penalty term.")
        .add_value(non_symmetric, "non-symmetric", "non-symmetric weighted interior penalty DG method")
        .add_value(incomplete,    "incomplete",    "incomplete weighted interior penalty DG method")
        .add_value(symmetric,     "symmetric",     "symmetric weighted interior penalty DG method")
        .close();
}

/*
*  Should be removed
template<class Model>
const Selection & TransportDG<Model>::EqData::get_output_selection() {
    return Model::ModelEqData::get_output_selection_input_type(
            "DG",
            "Implicit in time Discontinuous Galerkin solver")
        .copy_values(EqData().make_output_field_selection("").close())
        ConvectionTransport::EqData().output_fields
                                    .make_output_field_selection(
                                        "ConvectionTransport_output_fields",
                                        "Selection of output fields for Convection Solute Transport model.")
                                    .close()),
        .close();
}
*/

template<class Model>
const Record & TransportDG<Model>::get_input_type() {
    std::string equation_name = std::string(Model::ModelEqData::name()) + "_DG";
    return Model::get_input_type("DG", "Discontinuous Galerkin (DG) solver")
        .declare_key("solver", LinSys_PETSC::get_input_type(), Default("{}"),
                "Solver for the linear system.")
        .declare_key("user_fields", Array(
                TransportDG<Model>::EqFields()
                    .make_user_field_type(equation_name)),
                IT::Default::optional(),
                "Input fields of the equation defined by user.")
        .declare_key("input_fields", Array(
                TransportDG<Model>::EqFields()
                    .make_field_descriptor_type(equation_name)),
                IT::Default::obligatory(),
                "Input fields of the equation.")
        .declare_key("dg_variant", TransportDG<Model>::get_dg_variant_selection_input_type(), Default("\"non-symmetric\""),
                "Variant of the interior penalty discontinuous Galerkin method.")
        .declare_key("dg_order", Integer(0,3), Default("1"),
                "Polynomial order for the finite element in DG method (order 0 is suitable if there is no diffusion/dispersion).")
        .declare_key("init_projection", Bool(), Default("true"),
                "If true, use DG projection of the initial condition field."
                "Otherwise, evaluate initial condition field directly (well suited for reading native data).")
        .declare_key("output",
                EqFields().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ " + Model::ModelEqData::default_output_field() + "] }"),
                "Specification of output fields and output times.")
        .close();
}

template<class Model>
const int TransportDG<Model>::registrar =
        Input::register_class< TransportDG<Model>, Mesh &, const Input::Record>(std::string(Model::ModelEqData::name()) + "_DG") +
        TransportDG<Model>::get_input_type().size();



template<class Model>
TransportDG<Model>::EqFields::EqFields() : Model::ModelEqFields()
{
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of diffusive transfer through fractures (for each substance).")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(FieldFlag::in_main_matrix);

    *this+=dg_penalty
            .name("dg_penalty")
            .description(
            "Penalty parameter influencing the discontinuity of the solution (for each substance). "
            "Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);

    *this += region_id.name("region_id")
                .units( UnitSI::dimensionless())
                .flags(FieldFlag::equation_external_output)
                .description("Region ids.");
                
    *this += subdomain.name("subdomain")
      .units( UnitSI::dimensionless() )
      .flags(FieldFlag::equation_external_output)
      .description("Subdomain ids of the domain decomposition.");


    // add all input fields to the output list
    output_fields += *this;

    this->add_coords_field();
    this->set_default_fieldset();
}



// return the ratio of longest and shortest edge
template<class Model>
double TransportDG<Model>::EqData::elem_anisotropy(ElementAccessor<3> e) const
{
    double h_max = 0, h_min = numeric_limits<double>::infinity();
    for (unsigned int i=0; i<e->n_nodes(); i++)
        for (unsigned int j=i+1; j<e->n_nodes(); j++)
        {
            double dist = arma::norm(*e.node(i) - *e.node(j));
            h_max = max(h_max, dist);
            h_min = min(h_min, dist);
        }
    return h_max/h_min;
}



template<class Model>
void TransportDG<Model>::EqData::set_DG_parameters_boundary(Side side,
            const int K_size,
            const vector<arma::mat33> &K,
            const double flux,
            const arma::vec3 &normal_vector,
            const double alpha,
            double &gamma)
{
    double delta = 0, h = 0;

    // calculate the side diameter
    if (side.dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (unsigned int i=0; i<side.n_nodes(); i++)
            for (unsigned int j=i+1; j<side.n_nodes(); j++) {
                double dist = arma::norm(*side.node(i) - *side.node(j));
                h = max(h, dist);
            }

    }

    // delta is set to the average value of Kn.n on the side
    for (int k=0; k<K_size; k++)
        delta += dot(K[k]*normal_vector,normal_vector);
    delta /= K_size;

    gamma = 0.5*fabs(flux) + alpha/h*delta*elem_anisotropy(side.element());
}



template<typename Model>
TransportDG<Model>::TransportDG(Mesh & init_mesh, const Input::Record in_rec)
        : Model(init_mesh, in_rec),
          input_rec(in_rec),
          allocation_done(false),
          mass_assembly_(nullptr)
{
    // Can not use name() + "constructor" here, since START_TIMER only accepts const char *
    // due to constexpr optimization.
    START_TIMER(Model::ModelEqData::name());
    // Check that Model is derived from AdvectionDiffusionModel.
    static_assert(std::is_base_of<AdvectionDiffusionModel, Model>::value, "");

    eq_data_ = make_shared<EqData>();
    eq_fields_ = make_shared<EqFields>();
    this->eq_fieldset_ = eq_fields_.get();
    Model::init_balance(in_rec);


    // Set up physical parameters.
    eq_fields_->set_mesh(init_mesh);
    eq_fields_->region_id = GenericField<3>::region_id(*Model::mesh_);
    eq_fields_->subdomain = GenericField<3>::subdomain(*Model::mesh_);


    // DG data parameters
    eq_data_->dg_variant = in_rec.val<DGVariant>("dg_variant");
    eq_data_->dg_order = in_rec.val<unsigned int>("dg_order");
    
    Model::init_from_input(in_rec);

	MixedPtr<FE_P_disc> fe(eq_data_->dg_order);
	shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(Model::mesh_, fe);
	eq_data_->dh_ = make_shared<DOFHandlerMultiDim>(*Model::mesh_);
	eq_data_->dh_->distribute_dofs(ds);
    //DebugOut().fmt("TDG: solution size {}\n", eq_data_->dh_->n_global_dofs());

}


template<class Model>
void TransportDG<Model>::initialize()
{
    eq_fields_->set_components(eq_data_->substances_.names());
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), *(Model::time_) );
    eq_data_->set_time_governor(Model::time_);
    eq_data_->balance_ = this->balance();
    eq_fields_->initialize();

    // DG stabilization parameters on boundary edges
    eq_data_->gamma.resize(eq_data_->n_substances());
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        eq_data_->gamma[sbi].resize(Model::mesh_->boundary_.size());

    // Resize coefficient arrays
    eq_data_->max_edg_sides = max(Model::mesh_->max_edge_sides(1), max(Model::mesh_->max_edge_sides(2), Model::mesh_->max_edge_sides(3)));
    ret_sources.resize(eq_data_->n_substances());
    ret_sources_prev.resize(eq_data_->n_substances());

    Input::Array user_fields_arr;
    if (input_rec.opt_val("user_fields", user_fields_arr)) {
       	eq_fields_->init_user_fields(user_fields_arr, Model::time_->step().end(), eq_fields_->output_fields);
    }

    eq_data_->output_vec.resize(eq_data_->n_substances());
    eq_fields_->output_field.set_components(eq_data_->substances_.names());
    eq_fields_->output_field.set_mesh(*Model::mesh_);
    eq_fields_->output_fields.set_mesh(*Model::mesh_);
    eq_fields_->output_type(OutputTime::CORNER_DATA);

    eq_fields_->output_field.setup_components();
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
        // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
        auto output_field_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_);
        eq_fields_->output_field[sbi].set(output_field_ptr, 0);
        eq_data_->output_vec[sbi] = output_field_ptr->vec();
    }

    // set time marks for writing the output
    eq_fields_->output_fields.initialize(Model::output_stream_, Model::mesh_, input_rec.val<Input::Record>("output"), this->time());

    // equation default PETSc solver options
    std::string petsc_default_opts;
    if (eq_data_->dh_->distr()->np() == 1)
      petsc_default_opts = "-ksp_type bcgs -pc_type ilu -pc_factor_levels 2 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    else
      petsc_default_opts = "-ksp_type bcgs -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3 -sub_pc_factor_fill 6.0";
    
    // allocate matrix and vector structures
    eq_data_->ls    = new LinSys*[eq_data_->n_substances()];
    eq_data_->ls_dt = new LinSys*[eq_data_->n_substances()];
    eq_data_->conc_fe.resize(eq_data_->n_substances());
    
    MixedPtr<FE_P_disc> fe(0);
    shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(Model::mesh_, fe);
    eq_data_->dh_p0 = make_shared<DOFHandlerMultiDim>(*Model::mesh_);
    eq_data_->dh_p0->distribute_dofs(ds);

    stiffness_matrix.resize(eq_data_->n_substances(), nullptr);
    mass_matrix.resize(eq_data_->n_substances(), nullptr);
    rhs.resize(eq_data_->n_substances(), nullptr);
    mass_vec.resize(eq_data_->n_substances(), nullptr);
    eq_data_->ret_vec.resize(eq_data_->n_substances(), nullptr);

    for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++) {
        eq_data_->ls[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
        eq_data_->ls[sbi]->set_solution(eq_data_->output_vec[sbi].petsc_vec());

        eq_data_->ls_dt[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)eq_data_->ls_dt[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
        
        eq_data_->conc_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_p0);
        
        VecDuplicate(eq_data_->ls[sbi]->get_solution(), &eq_data_->ret_vec[sbi]);
    }


    init_projection = input_rec.val<bool>("init_projection");

    // create assemblation object, finite element structures and distribute DOFs
	mass_assembly_ = new GenericAssembly< MassAssemblyDim >(eq_fields_.get(), eq_data_.get());
	stiffness_assembly_ = new GenericAssembly< StiffnessAssemblyDim >(eq_fields_.get(), eq_data_.get());
	sources_assembly_ = new GenericAssembly< SourcesAssemblyDim >(eq_fields_.get(), eq_data_.get());
	bdr_cond_assembly_ = new GenericAssembly< BdrConditionAssemblyDim >(eq_fields_.get(), eq_data_.get());
    
    if(init_projection)
	    init_assembly_ = new GenericAssembly< InitProjectionAssemblyDim >(eq_fields_.get(), eq_data_.get());
    else
        init_assembly_ = new GenericAssembly< InitConditionAssemblyDim >(eq_fields_.get(), eq_data_.get());

    // initialization of balance object
    Model::balance_->allocate(eq_data_->dh_->distr()->lsize(), mass_assembly_->eval_points()->max_size());

    int qsize = mass_assembly_->eval_points()->max_size();
    eq_data_->dif_coef.resize(eq_data_->n_substances());
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
        eq_data_->dif_coef[sbi].resize(qsize);
    }

    eq_fields_->init_condition.setup_components();
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
    	eq_fields_->init_condition[sbi].add_factory( std::make_shared<FieldFE<3, FieldValue<3>::Scalar>::NativeFactory>(sbi, eq_data_->dh_));
    }
}


template<class Model>
TransportDG<Model>::~TransportDG()
{
    delete Model::time_;

    if (eq_data_->gamma.size() > 0) {
        // initialize called

        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
            if (eq_data_->ls != nullptr) {
                delete eq_data_->ls[i];
                delete eq_data_->ls_dt[i];
            }

            if (stiffness_matrix.size() > 0) {
                if (stiffness_matrix[i])
                    chkerr(MatDestroy(&stiffness_matrix[i]));
                if (mass_matrix[i])
                    chkerr(MatDestroy(&mass_matrix[i]));
                if (rhs[i])
                	chkerr(VecDestroy(&rhs[i]));
                if (mass_vec[i])
                	chkerr(VecDestroy(&mass_vec[i]));
                if (eq_data_->ret_vec[i])
                	chkerr(VecDestroy(&eq_data_->ret_vec[i]));
            }
        }
        if (eq_data_->ls != nullptr) {
            delete[] eq_data_->ls;
            delete[] eq_data_->ls_dt;
            eq_data_->ls = nullptr;
        }
        //delete[] stiffness_matrix;
        //delete[] mass_matrix;
        //delete[] rhs;
        //delete[] mass_vec;
        //delete[] ret_vec;

        if (mass_assembly_ != nullptr) {
            delete mass_assembly_;
            delete stiffness_assembly_;
            delete sources_assembly_;
            delete bdr_cond_assembly_;
            delete init_assembly_;
        }
    }


}


template<class Model>
void TransportDG<Model>::zero_time_step()
{
    START_TIMER(Model::ModelEqData::name());
    eq_fields_->mark_input_times( *(Model::time_) );
    eq_fields_->set_time(Model::time_->step(), LimitSide::left);
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "transport DG") ) {
        WarningOut() << ss.str();
    }


    // set initial conditions
    set_initial_condition();
    for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++)
        ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_initial_guess_nonzero();

    // check first time assembly - needs preallocation
    if (!allocation_done) preallocate();

    // after preallocation we assemble the matrices and vectors required for mass balance
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
    {
        Model::balance_->calculate_instant(eq_data_->subst_idx_[sbi], eq_data_->ls[sbi]->get_solution());

        // add sources due to sorption
        ret_sources_prev[sbi] = 0;
    }

    output_data();
}


template<class Model>
void TransportDG<Model>::preallocate()
{
    // preallocate system matrix
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        // preallocate system matrix
    	eq_data_->ls[i]->start_allocation();
        stiffness_matrix[i] = NULL;
        rhs[i] = NULL;

        // preallocate mass matrix
        eq_data_->ls_dt[i]->start_allocation();
        mass_matrix[i] = NULL;
        VecZeroEntries(eq_data_->ret_vec[i]);
    }
    stiffness_assembly_->assemble(eq_data_->dh_);
    mass_assembly_->assemble(eq_data_->dh_);
    sources_assembly_->assemble(eq_data_->dh_);
    bdr_cond_assembly_->assemble(eq_data_->dh_);
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
      VecAssemblyBegin(eq_data_->ret_vec[i]);
      VecAssemblyEnd(eq_data_->ret_vec[i]);
    }

    allocation_done = true;
}



template<class Model>
void TransportDG<Model>::update_solution()
{
    START_TIMER("DG-ONE STEP");

    Model::time_->next_time();
    Model::time_->view("TDG");
    
    START_TIMER("data reinit");
    eq_fields_->set_time(Model::time_->step(), LimitSide::left);
    END_TIMER("data reinit");

    // assemble mass matrix
    if (mass_matrix[0] == NULL || eq_fields_->subset(FieldFlag::in_time_term).changed() )
    {
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
        	eq_data_->ls_dt[i]->start_add_assembly();
        	eq_data_->ls_dt[i]->mat_zero_entries();
            VecZeroEntries(eq_data_->ret_vec[i]);
        }
        mass_assembly_->assemble(eq_data_->dh_);
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
        	eq_data_->ls_dt[i]->finish_assembly();
            VecAssemblyBegin(eq_data_->ret_vec[i]);
            VecAssemblyEnd(eq_data_->ret_vec[i]);
            // construct mass_vec for initial time
            if (mass_matrix[i] == NULL)
            {
                VecDuplicate(eq_data_->ls[i]->get_solution(), &mass_vec[i]);
                MatMult(*(eq_data_->ls_dt[i]->get_matrix()), eq_data_->ls[i]->get_solution(), mass_vec[i]);
                MatConvert(*( eq_data_->ls_dt[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &mass_matrix[i]);
            }
            else
                MatCopy(*( eq_data_->ls_dt[i]->get_matrix() ), mass_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble stiffness matrix
    if (stiffness_matrix[0] == NULL
            || eq_fields_->subset(FieldFlag::in_main_matrix).changed()
            || eq_fields_->flow_flux.changed())
    {
        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
            eq_data_->ls[i]->start_add_assembly();
            eq_data_->ls[i]->mat_zero_entries();
        }
        stiffness_assembly_->assemble(eq_data_->dh_);
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
        	eq_data_->ls[i]->finish_assembly();

            if (stiffness_matrix[i] == NULL)
                MatConvert(*( eq_data_->ls[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
            else
                MatCopy(*( eq_data_->ls[i]->get_matrix() ), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (rhs[0] == NULL
            || eq_fields_->subset(FieldFlag::in_rhs).changed()
            || eq_fields_->flow_flux.changed())
    {
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
            eq_data_->ls[i]->start_add_assembly();
            eq_data_->ls[i]->rhs_zero_entries();
        }
        sources_assembly_->assemble(eq_data_->dh_);
        bdr_cond_assembly_->assemble(eq_data_->dh_);
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
            eq_data_->ls[i]->finish_assembly();

            if (rhs[i] == nullptr) VecDuplicate(*( eq_data_->ls[i]->get_rhs() ), &rhs[i]);
            VecCopy(*( eq_data_->ls[i]->get_rhs() ), rhs[i]);
        }
    }

    /* Apply backward Euler time integration.
    *
    * Denoting A the stiffness matrix and M the mass matrix, the algebraic system at the k-th time level reads
    *
    *   (1/dt M + A)u^k = f + 1/dt M.u^{k-1}
    *
    * Hence we modify at each time level the right hand side:
    *
    *   f^k = f + 1/dt M u^{k-1},
    *
    * where f stands for the term stemming from the force and boundary conditions.
    * Accordingly, we set
    *
    *   A^k = A + 1/dt M.
    *
    */
    Mat m;
    START_TIMER("solve");
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        MatConvert(stiffness_matrix[i], MATSAME, MAT_INITIAL_MATRIX, &m);
        MatAXPY(m, 1./Model::time_->dt(), mass_matrix[i], SUBSET_NONZERO_PATTERN);
        eq_data_->ls[i]->set_matrix(m, DIFFERENT_NONZERO_PATTERN);
        Vec w;
        VecDuplicate(rhs[i], &w);
        VecWAXPY(w, 1./Model::time_->dt(), mass_vec[i], rhs[i]);
        eq_data_->ls[i]->set_rhs(w);

        VecDestroy(&w);
        chkerr(MatDestroy(&m));

        eq_data_->ls[i]->solve();

        // update mass_vec due to possible changes in mass matrix
        MatMult(*(eq_data_->ls_dt[i]->get_matrix()), eq_data_->ls[i]->get_solution(), mass_vec[i]);
    }
    END_TIMER("solve");

    // Possibly output matrices for debug reasons.
    // for (unsigned int i=0; i<eq_data_->n_substances(); i++){
    //     string conc_name = eq_data_->substances().names()[i] + "_" + std::to_string(eq_data_->time_->step().index());
    //     eq_data_->ls[i]->view("stiff_" + conc_name);
    //     eq_data_->ls_dt[i]->view("mass_" + conc_name);
    // }

    calculate_cumulative_balance();

    END_TIMER("DG-ONE STEP");
}


template<class Model>
void TransportDG<Model>::compute_p0_interpolation()
{
    // calculate element averages of solution
	for (auto cell : eq_data_->dh_->own_range() )
    {
		LocDofVec loc_dof_indices = cell.get_loc_dof_indices();
		unsigned int n_dofs=loc_dof_indices.n_rows;

        DHCellAccessor dh_p0_cell = eq_data_->dh_p0->cell_accessor_from_element(cell.elm_idx());
        IntIdx dof_p0 = dh_p0_cell.get_loc_dof_indices()[0];

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
        {
            eq_data_->conc_fe[sbi]->vec().set(dof_p0, 0);

            for (unsigned int j=0; j<n_dofs; ++j)
                eq_data_->conc_fe[sbi]->vec().add( dof_p0, eq_data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]] );

            eq_data_->conc_fe[sbi]->vec().set( dof_p0, max(eq_data_->conc_fe[sbi]->vec().get(dof_p0)/n_dofs, 0.) );
        }
    }
}




template<class Model>
void TransportDG<Model>::output_data()
{
    //if (!Model::time_->is_current( Model::time_->marks().type_output() )) return;


    START_TIMER("DG-OUTPUT");

    // gather the solution from all processors
    eq_fields_->output_fields.set_time( this->time().step(), LimitSide::left);
    //if (eq_fields_->output_fields.is_field_output_time(eq_fields_->output_field, this->time().step()) )
    eq_fields_->output_fields.output(this->time().step());

    Model::output_data();
    
    START_TIMER("TOS-balance");
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
      Model::balance_->calculate_instant(eq_data_->subst_idx_[sbi], eq_data_->ls[sbi]->get_solution());
    Model::balance_->output();
    END_TIMER("TOS-balance");

    END_TIMER("DG-OUTPUT");
}


template<class Model>
void TransportDG<Model>::calculate_cumulative_balance()
{
    if (Model::balance_->cumulative())
    {
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
        {
            Model::balance_->calculate_cumulative(eq_data_->subst_idx_[sbi], eq_data_->ls[sbi]->get_solution());

            // update source increment due to retardation
            VecDot(eq_data_->ret_vec[sbi], eq_data_->ls[sbi]->get_solution(), &ret_sources[sbi]);

            Model::balance_->add_cumulative_source(eq_data_->subst_idx_[sbi], (ret_sources[sbi]-ret_sources_prev[sbi])/Model::time_->dt());
            ret_sources_prev[sbi] = ret_sources[sbi];
        }
    }
}




template<class Model>
void TransportDG<Model>::set_initial_condition()
{
    if(init_projection)
    {
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
            eq_data_->ls[sbi]->start_allocation();
        
        init_assembly_->assemble(eq_data_->dh_);

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
            eq_data_->ls[sbi]->start_add_assembly();

        init_assembly_->assemble(eq_data_->dh_);

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        {
            eq_data_->ls[sbi]->finish_assembly();
            eq_data_->ls[sbi]->solve();
        }
    }
    else
        init_assembly_->assemble(eq_data_->dh_);
}


template<class Model>
void TransportDG<Model>::update_after_reactions(bool solution_changed)
{
    if (solution_changed)
    {
    	for (auto cell : eq_data_->dh_->own_range() )
        {
    	    LocDofVec loc_dof_indices = cell.get_loc_dof_indices();
            unsigned int n_dofs=loc_dof_indices.n_rows;
            
            DHCellAccessor dh_p0_cell = eq_data_->dh_p0->cell_accessor_from_element(cell.elm_idx());
            IntIdx dof_p0 = dh_p0_cell.get_loc_dof_indices()[0];

            for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
            {
                double old_average = 0;
                for (unsigned int j=0; j<n_dofs; ++j)
                    old_average += eq_data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]];
                old_average /= n_dofs;

                for (unsigned int j=0; j<n_dofs; ++j)
                    eq_data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]]
                        += eq_data_->conc_fe[sbi]->vec().get(dof_p0) - old_average;
            }
        }
    }
    // update mass_vec for the case that mass matrix changes in next time step
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
        MatMult(*(eq_data_->ls_dt[sbi]->get_matrix()), eq_data_->ls[sbi]->get_solution(), mass_vec[sbi]);
}






template class TransportDG<ConcentrationTransportModel>;
template class TransportDG<HeatTransferModel>;




