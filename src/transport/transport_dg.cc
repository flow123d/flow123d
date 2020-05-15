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
        .declare_key("input_fields", Array(
                TransportDG<Model>::EqData()
                    .make_field_descriptor_type(equation_name)),
                IT::Default::obligatory(),
                "Input fields of the equation.")
        .declare_key("dg_variant", TransportDG<Model>::get_dg_variant_selection_input_type(), Default("\"non-symmetric\""),
                "Variant of the interior penalty discontinuous Galerkin method.")
        .declare_key("dg_order", Integer(0,3), Default("1"),
                "Polynomial order for the finite element in DG method (order 0 is suitable if there is no diffusion/dispersion).")
        .declare_key("output",
                EqData().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ " + Model::ModelEqData::default_output_field() + "] }"),
                "Specification of output fields and output times.")
        .close();
}

template<class Model>
const int TransportDG<Model>::registrar =
        Input::register_class< TransportDG<Model>, Mesh &, const Input::Record>(std::string(Model::ModelEqData::name()) + "_DG") +
        TransportDG<Model>::get_input_type().size();



template<class Model>
TransportDG<Model>::EqData::EqData() : Model::ModelEqData()
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
          allocation_done(false)
{
    // Can not use name() + "constructor" here, since START_TIMER only accepts const char *
    // due to constexpr optimization.
    START_TIMER(Model::ModelEqData::name());
    // Check that Model is derived from AdvectionDiffusionModel.
    static_assert(std::is_base_of<AdvectionDiffusionModel, Model>::value, "");

    data_ = make_shared<EqData>();
    this->eq_data_ = data_.get();


    // Set up physical parameters.
    data_->set_mesh(init_mesh);
    data_->region_id = GenericField<3>::region_id(*Model::mesh_);
    data_->subdomain = GenericField<3>::subdomain(*Model::mesh_);


    // DG variant and order
    data_->dg_variant = in_rec.val<DGVariant>("dg_variant");
    data_->dg_order = in_rec.val<unsigned int>("dg_order");
    
    Model::init_from_input(in_rec);

    // create assemblation object, finite element structures and distribute DOFs
	data_->mass_assembly_ = new GenericAssembly< MassAssemblyDim >(data_.get(), ActiveIntegrals::bulk);
	data_->stiffness_assembly_ = new GenericAssembly< StiffnessAssemblyDim >(data_.get(),
	        (ActiveIntegrals::bulk | ActiveIntegrals::edge | ActiveIntegrals::coupling | ActiveIntegrals::boundary) );
	data_->sources_assembly_ = new GenericAssembly< SourcesAssemblyDim >(data_.get(), ActiveIntegrals::bulk);
	data_->bdr_cond_assembly_ = new GenericAssembly< BdrConditionAssemblyDim >(data_.get(), ActiveIntegrals::bulk);
	data_->init_cond_assembly_ = new GenericAssembly< InitConditionAssemblyDim >(data_.get(), ActiveIntegrals::bulk);

	MixedPtr<FE_P_disc> fe(data_->dg_order);
	shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(Model::mesh_, fe);
	data_->dh_ = make_shared<DOFHandlerMultiDim>(*Model::mesh_);
	data_->dh_->distribute_dofs(ds);
    //DebugOut().fmt("TDG: solution size {}\n", data_->dh_->n_global_dofs());

}


template<class Model>
void TransportDG<Model>::initialize()
{
    data_->set_components(Model::substances_.names());
    data_->set_input_list( input_rec.val<Input::Array>("input_fields"), *(Model::time_) );

    // DG stabilization parameters on boundary edges
    data_->gamma.resize(Model::n_substances());
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    	data_->gamma[sbi].resize(Model::mesh_->boundary_.size());

    // Resize coefficient arrays
    int qsize = data_->mass_assembly_->eval_points()->max_size();
    int max_edg_sides = max(Model::mesh_->max_edge_sides(1), max(Model::mesh_->max_edge_sides(2), Model::mesh_->max_edge_sides(3)));
    ret_sources.resize(Model::n_substances());
    ret_sources_prev.resize(Model::n_substances());
    data_->ad_coef.resize(Model::n_substances());
    data_->dif_coef.resize(Model::n_substances());
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    {
        data_->ad_coef[sbi].resize(qsize);
        data_->dif_coef[sbi].resize(qsize);
    }
    data_->ad_coef_edg.resize(max_edg_sides);
    data_->dif_coef_edg.resize(max_edg_sides);
    for (int sd=0; sd<max_edg_sides; sd++)
    {
        data_->ad_coef_edg[sd].resize(Model::n_substances());
        data_->dif_coef_edg[sd].resize(Model::n_substances());
        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
            data_->ad_coef_edg[sd][sbi].resize(qsize);
            data_->dif_coef_edg[sd][sbi].resize(qsize);
        }
    }

    output_vec.resize(Model::n_substances());
    data_->output_field.set_components(Model::substances_.names());
    data_->output_field.set_mesh(*Model::mesh_);
    data_->output_type(OutputTime::CORNER_DATA);

    data_->output_field.setup_components();
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    {
        // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_field_ptr(new FieldFE<3, FieldValue<3>::Scalar>);
        output_vec[sbi] = output_field_ptr->set_fe_data(data_->dh_);
        data_->output_field[sbi].set_field(Model::mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
    }

    // set time marks for writing the output
    data_->output_fields.initialize(Model::output_stream_, Model::mesh_, input_rec.val<Input::Record>("output"), this->time());

    // equation default PETSc solver options
    std::string petsc_default_opts;
    if (data_->dh_->distr()->np() == 1)
      petsc_default_opts = "-ksp_type bcgs -pc_type ilu -pc_factor_levels 2 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    else
      petsc_default_opts = "-ksp_type bcgs -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3 -sub_pc_factor_fill 6.0";
    
    // allocate matrix and vector structures
    data_->ls    = new LinSys*[Model::n_substances()];
    data_->ls_dt = new LinSys*[Model::n_substances()];
    solution_elem_ = new double*[Model::n_substances()];

    stiffness_matrix.resize(Model::n_substances(), nullptr);
    mass_matrix.resize(Model::n_substances(), nullptr);
    rhs.resize(Model::n_substances(), nullptr);
    mass_vec.resize(Model::n_substances(), nullptr);
    data_->ret_vec.resize(Model::n_substances(), nullptr);

    for (unsigned int sbi = 0; sbi < Model::n_substances(); sbi++) {
    	data_->ls[sbi] = new LinSys_PETSC(data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)data_->ls[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
        data_->ls[sbi]->set_solution(output_vec[sbi].petsc_vec());

        data_->ls_dt[sbi] = new LinSys_PETSC(data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)data_->ls_dt[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
        solution_elem_[sbi] = new double[Model::mesh_->get_el_ds()->lsize()];
        
        VecDuplicate(data_->ls[sbi]->get_solution(), &data_->ret_vec[sbi]);
    }


    // initialization of balance object
    Model::balance_->allocate(data_->dh_->distr()->lsize(), data_->mass_assembly_->eval_points()->max_size());

    // initialization of assembly object
    data_->mass_assembly_->multidim_assembly()[1_d]->initialize(*this);
    data_->mass_assembly_->multidim_assembly()[2_d]->initialize(*this);
    data_->mass_assembly_->multidim_assembly()[3_d]->initialize(*this);
    data_->stiffness_assembly_->multidim_assembly()[1_d]->initialize(*this);
    data_->stiffness_assembly_->multidim_assembly()[2_d]->initialize(*this);
    data_->stiffness_assembly_->multidim_assembly()[3_d]->initialize(*this);
    data_->sources_assembly_->multidim_assembly()[1_d]->initialize(*this);
    data_->sources_assembly_->multidim_assembly()[2_d]->initialize(*this);
    data_->sources_assembly_->multidim_assembly()[3_d]->initialize(*this);
    data_->bdr_cond_assembly_->multidim_assembly()[1_d]->initialize(*this);
    data_->bdr_cond_assembly_->multidim_assembly()[2_d]->initialize(*this);
    data_->bdr_cond_assembly_->multidim_assembly()[3_d]->initialize(*this);
    data_->init_cond_assembly_->multidim_assembly()[1_d]->initialize(*this);
    data_->init_cond_assembly_->multidim_assembly()[2_d]->initialize(*this);
    data_->init_cond_assembly_->multidim_assembly()[3_d]->initialize(*this);
}


template<class Model>
TransportDG<Model>::~TransportDG()
{
    delete Model::time_;

    if (data_->gamma.size() > 0) {
        // initialize called

        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            delete data_->ls[i];
            delete[] solution_elem_[i];
            delete data_->ls_dt[i];

            if (stiffness_matrix[i])
                chkerr(MatDestroy(&stiffness_matrix[i]));
            if (mass_matrix[i])
                chkerr(MatDestroy(&mass_matrix[i]));
            if (rhs[i])
            	chkerr(VecDestroy(&rhs[i]));
            if (mass_vec[i])
            	chkerr(VecDestroy(&mass_vec[i]));
            if (data_->ret_vec[i])
            	chkerr(VecDestroy(&data_->ret_vec[i]));
        }
        delete[] data_->ls;
        delete[] solution_elem_;
        delete[] data_->ls_dt;
        //delete[] stiffness_matrix;
        //delete[] mass_matrix;
        //delete[] rhs;
        //delete[] mass_vec;
        //delete[] ret_vec;

        delete data_->mass_assembly_;
        delete data_->stiffness_assembly_;
        delete data_->sources_assembly_;
        delete data_->bdr_cond_assembly_;
        delete data_->init_cond_assembly_;
    }

}


template<class Model>
void TransportDG<Model>::zero_time_step()
{
    START_TIMER(Model::ModelEqData::name());
    data_->mark_input_times( *(Model::time_) );
    data_->set_time(Model::time_->step(), LimitSide::left);
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "transport DG") ) {
        WarningOut() << ss.str();
    }


    // set initial conditions
    set_initial_condition();
    for (unsigned int sbi = 0; sbi < Model::n_substances(); sbi++)
        ( (LinSys_PETSC *)data_->ls[sbi] )->set_initial_guess_nonzero();

    // check first time assembly - needs preallocation
    if (!allocation_done) preallocate();

    // after preallocation we assemble the matrices and vectors required for mass balance
    for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
    {
        Model::balance_->calculate_instant(Model::subst_idx[sbi], data_->ls[sbi]->get_solution());

        // add sources due to sorption
        ret_sources_prev[sbi] = 0;
    }

    output_data();
}


template<class Model>
void TransportDG<Model>::preallocate()
{
    // preallocate system matrix
    for (unsigned int i=0; i<Model::n_substances(); i++)
    {
        // preallocate system matrix
    	data_->ls[i]->start_allocation();
        stiffness_matrix[i] = NULL;
        rhs[i] = NULL;

        // preallocate mass matrix
        data_->ls_dt[i]->start_allocation();
        mass_matrix[i] = NULL;
        VecZeroEntries(data_->ret_vec[i]);
    }
    START_TIMER("assemble_stiffness");
    data_->stiffness_assembly_->assemble(data_->dh_);
    END_TIMER("assemble_stiffness");
    START_TIMER("assemble_mass");
    data_->mass_assembly_->assemble(data_->dh_);
    END_TIMER("assemble_mass");
    START_TIMER("assemble_sources");
    data_->sources_assembly_->assemble(data_->dh_);
    END_TIMER("assemble_sources");
    START_TIMER("assemble_bc");
    data_->bdr_cond_assembly_->assemble(data_->dh_);
    END_TIMER("assemble_bc");
    for (unsigned int i=0; i<Model::n_substances(); i++)
    {
      VecAssemblyBegin(data_->ret_vec[i]);
      VecAssemblyEnd(data_->ret_vec[i]);
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
    data_->set_time(Model::time_->step(), LimitSide::left);
    END_TIMER("data reinit");
    
    // assemble mass matrix
    if (mass_matrix[0] == NULL || data_->subset(FieldFlag::in_time_term).changed() )
    {
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
        	data_->ls_dt[i]->start_add_assembly();
        	data_->ls_dt[i]->mat_zero_entries();
            VecZeroEntries(data_->ret_vec[i]);
        }
        START_TIMER("assemble_mass");
        data_->mass_assembly_->assemble(data_->dh_);
        END_TIMER("assemble_mass");
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
        	data_->ls_dt[i]->finish_assembly();
            VecAssemblyBegin(data_->ret_vec[i]);
            VecAssemblyEnd(data_->ret_vec[i]);
            // construct mass_vec for initial time
            if (mass_matrix[i] == NULL)
            {
                VecDuplicate(data_->ls[i]->get_solution(), &mass_vec[i]);
                MatMult(*(data_->ls_dt[i]->get_matrix()), data_->ls[i]->get_solution(), mass_vec[i]);
                MatConvert(*( data_->ls_dt[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &mass_matrix[i]);
            }
            else
                MatCopy(*( data_->ls_dt[i]->get_matrix() ), mass_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble stiffness matrix
    if (stiffness_matrix[0] == NULL
            || data_->subset(FieldFlag::in_main_matrix).changed()
            || Model::flux_changed)
    {
        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            data_->ls[i]->start_add_assembly();
            data_->ls[i]->mat_zero_entries();
        }
        START_TIMER("assemble_stiffness");
        data_->stiffness_assembly_->assemble(data_->dh_);
        END_TIMER("assemble_stiffness");
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
        	data_->ls[i]->finish_assembly();

            if (stiffness_matrix[i] == NULL)
                MatConvert(*( data_->ls[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
            else
                MatCopy(*( data_->ls[i]->get_matrix() ), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (rhs[0] == NULL
            || data_->subset(FieldFlag::in_rhs).changed()
            || Model::flux_changed)
    {
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            data_->ls[i]->start_add_assembly();
            data_->ls[i]->rhs_zero_entries();
        }
        START_TIMER("assemble_sources");
        data_->sources_assembly_->assemble(data_->dh_);
        END_TIMER("assemble_sources");
        START_TIMER("assemble_bc");
        data_->bdr_cond_assembly_->assemble(data_->dh_);
        END_TIMER("assemble_bc");
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            data_->ls[i]->finish_assembly();

            if (rhs[i] == nullptr) VecDuplicate(*( data_->ls[i]->get_rhs() ), &rhs[i]);
            VecCopy(*( data_->ls[i]->get_rhs() ), rhs[i]);
        }
    }
    
    Model::flux_changed = false;


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
    for (unsigned int i=0; i<Model::n_substances(); i++)
    {
        MatConvert(stiffness_matrix[i], MATSAME, MAT_INITIAL_MATRIX, &m);
        MatAXPY(m, 1./Model::time_->dt(), mass_matrix[i], SUBSET_NONZERO_PATTERN);
        data_->ls[i]->set_matrix(m, DIFFERENT_NONZERO_PATTERN);
        Vec w;
        VecDuplicate(rhs[i], &w);
        VecWAXPY(w, 1./Model::time_->dt(), mass_vec[i], rhs[i]);
        data_->ls[i]->set_rhs(w);

        VecDestroy(&w);
        chkerr(MatDestroy(&m));

        data_->ls[i]->solve();

        // update mass_vec due to possible changes in mass matrix
        MatMult(*(data_->ls_dt[i]->get_matrix()), data_->ls[i]->get_solution(), mass_vec[i]);
    }
    END_TIMER("solve");

    calculate_cumulative_balance();

    END_TIMER("DG-ONE STEP");
}


template<class Model>
void TransportDG<Model>::calculate_concentration_matrix()
{
    // calculate element averages of solution
	unsigned int i_cell=0;
	for (auto cell : data_->dh_->own_range() )
    {
		LocDofVec loc_dof_indices = cell.get_loc_dof_indices();
		unsigned int n_dofs=loc_dof_indices.n_rows;

        for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
        {
            solution_elem_[sbi][i_cell] = 0;

            for (unsigned int j=0; j<n_dofs; ++j)
                solution_elem_[sbi][i_cell] += data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]];

            solution_elem_[sbi][i_cell] = max(solution_elem_[sbi][i_cell]/n_dofs, 0.);
        }
        ++i_cell;
    }
}




template<class Model>
void TransportDG<Model>::output_data()
{
    //if (!Model::time_->is_current( Model::time_->marks().type_output() )) return;


    START_TIMER("DG-OUTPUT");

    // gather the solution from all processors
    data_->output_fields.set_time( this->time().step(), LimitSide::left);
    //if (data_->output_fields.is_field_output_time(data_->output_field, this->time().step()) )
    data_->output_fields.output(this->time().step());

    Model::output_data();
    
    START_TIMER("TOS-balance");
    for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
      Model::balance_->calculate_instant(Model::subst_idx[sbi], data_->ls[sbi]->get_solution());
    Model::balance_->output();
    END_TIMER("TOS-balance");

    END_TIMER("DG-OUTPUT");
}


template<class Model>
void TransportDG<Model>::calculate_cumulative_balance()
{
    if (Model::balance_->cumulative())
    {
        for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
        {
            Model::balance_->calculate_cumulative(Model::subst_idx[sbi], data_->ls[sbi]->get_solution());

            // update source increment due to retardation
            VecDot(data_->ret_vec[sbi], data_->ls[sbi]->get_solution(), &ret_sources[sbi]);

            Model::balance_->add_cumulative_source(Model::subst_idx[sbi], (ret_sources[sbi]-ret_sources_prev[sbi])/Model::time_->dt());
            ret_sources_prev[sbi] = ret_sources[sbi];
        }
    }
}




template<class Model>
void TransportDG<Model>::set_initial_condition()
{
    START_TIMER("set_init_cond");
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        data_->ls[sbi]->start_allocation();
    data_->init_cond_assembly_->assemble(data_->dh_);

    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        data_->ls[sbi]->start_add_assembly();
    data_->init_cond_assembly_->assemble(data_->dh_);

    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    {
        data_->ls[sbi]->finish_assembly();
        data_->ls[sbi]->solve();
    }
    END_TIMER("set_init_cond");
}


template<class Model>
void TransportDG<Model>::get_par_info(LongIdx * &el_4_loc, Distribution * &el_ds)
{
    el_4_loc = Model::mesh_->get_el_4_loc();
    el_ds = Model::mesh_->get_el_ds();
}


template<class Model>
void TransportDG<Model>::update_after_reactions(bool solution_changed)
{
    if (solution_changed)
    {
    	unsigned int i_cell=0;
    	for (auto cell : data_->dh_->own_range() )
        {
    	    LocDofVec loc_dof_indices = cell.get_loc_dof_indices();
            unsigned int n_dofs=loc_dof_indices.n_rows;

            for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
            {
                double old_average = 0;
                for (unsigned int j=0; j<n_dofs; ++j)
                    old_average += data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]];
                old_average /= n_dofs;

                for (unsigned int j=0; j<n_dofs; ++j)
                    data_->ls[sbi]->get_solution_array()[loc_dof_indices[j]] += solution_elem_[sbi][i_cell] - old_average;
            }
            ++i_cell;
        }
    }
    // update mass_vec for the case that mass matrix changes in next time step
    for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
        MatMult(*(data_->ls_dt[sbi]->get_matrix()), data_->ls[sbi]->get_solution(), mass_vec[sbi]);
}

template<class Model>
LongIdx *TransportDG<Model>::get_row_4_el()
{
    return Model::mesh_->get_row_4_el();
}






template class TransportDG<ConcentrationTransportModel>;
template class TransportDG<HeatTransferModel>;




