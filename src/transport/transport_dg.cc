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

#include "system/sys_profiler.hh"
#include "transport/transport_dg.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/dh_cell_accessor.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"
#include "la/linsys_PETSC.hh"
#include "flow/mh_dofhandler.hh"
#include "coupling/balance.hh"
#include "transport/advection_diffusion_model.hh"
#include "transport/concentration_model.hh"
#include "transport/heat_model.hh"
#include "transport/assembly_dg.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "input/factory.hh"
#include "fields/equation_output.hh"
#include "mesh/long_idx.hh"
#include "mesh/accessors.hh"

FLOW123D_FORCE_LINK_IN_CHILD(concentrationTransportModel);
FLOW123D_FORCE_LINK_IN_CHILD(heatModel);



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
            h_max = max(h_max, e.node(i)->distance(*e.node(j)));
            h_min = min(h_min, e.node(i)->distance(*e.node(j)));
        }
    return h_max/h_min;
}



template<class Model>
void TransportDG<Model>::EqData::set_DG_parameters_boundary(const Side *side,
            const int K_size,
            const vector<arma::mat33> &K,
            const double flux,
            const arma::vec3 &normal_vector,
            const double alpha,
            double &gamma)
{
    double delta = 0, h = 0;

    // calculate the side diameter
    if (side->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (unsigned int i=0; i<side->n_nodes(); i++)
            for (unsigned int j=i+1; j<side->n_nodes(); j++)
                h = max(h, side->node(i)->distance( *side->node(j).node() ));
    }

    // delta is set to the average value of Kn.n on the side
    for (int k=0; k<K_size; k++)
        delta += dot(K[k]*normal_vector,normal_vector);
    delta /= K_size;

    gamma = 0.5*fabs(flux) + alpha/h*delta*elem_anisotropy(side->element());
}





template<class Model>
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

    // create finite element structures and distribute DOFs
	assembly1_ = std::make_shared<AssemblyDG<1, Model>>(data_, *this);
	assembly2_ = std::make_shared<AssemblyDG<2, Model>>(data_, *this);
	assembly3_ = std::make_shared<AssemblyDG<3, Model>>(data_, *this);
	multidim_assembly_.push_back( std::dynamic_pointer_cast<AssemblyDGBase>(assembly1_) );
	multidim_assembly_.push_back( std::dynamic_pointer_cast<AssemblyDGBase>(assembly2_) );
	multidim_assembly_.push_back( std::dynamic_pointer_cast<AssemblyDGBase>(assembly3_) );
	shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(Model::mesh_, assembly1_->fe_low(), assembly1_->fe(), assembly2_->fe(), assembly3_->fe());
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
    int qsize = max(assembly1_->quad_low()->size(), max(assembly1_->quad()->size(), max(assembly2_->quad()->size(), assembly3_->quad()->size())));
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
    Model::balance_->allocate(data_->dh_->distr()->lsize(),
            max(assembly1_->fe()->n_dofs(), max(assembly2_->fe()->n_dofs(), assembly3_->fe()->n_dofs())));

    initialize_assembly_objects();
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
    assemble_stiffness_matrix();
    assemble_mass_matrix();
    set_sources();
    set_boundary_conditions();
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
        assemble_mass_matrix();
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
        assemble_stiffness_matrix();
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
        set_sources();
        set_boundary_conditions();
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

        unsigned int n_dofs;
        switch (cell.dim())
        {
        case 1:
            n_dofs = assembly1_->fe()->n_dofs();
            break;
        case 2:
            n_dofs = assembly2_->fe()->n_dofs();
            break;
        case 3:
            n_dofs = assembly3_->fe()->n_dofs();
            break;
        }

        std::vector<LongIdx> dof_indices(n_dofs);
        cell.get_dof_indices(dof_indices);

        for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
        {
            solution_elem_[sbi][i_cell] = 0;

            for (unsigned int j=0; j<n_dofs; ++j)
                solution_elem_[sbi][i_cell] += data_->ls[sbi]->get_solution_array()[dof_indices[j]-data_->dh_->distr()->begin()];

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
void TransportDG<Model>::assemble_mass_matrix()
{
    START_TIMER("assemble_mass");
    Model::balance_->start_mass_assembly(Model::subst_idx);
    for (auto cell : data_->dh_->own_range() )
    {
        multidim_assembly_[ cell.dim()-1 ]->assemble_mass_matrix(cell);
    }
    Model::balance_->finish_mass_assembly(Model::subst_idx);
    END_TIMER("assemble_mass");
}



template<class Model>
void TransportDG<Model>::assemble_stiffness_matrix()
{
  START_TIMER("assemble_stiffness");
    for (auto cell : data_->dh_->own_range() )
    {
        START_TIMER("assemble_volume_integrals");
        multidim_assembly_[ cell.dim()-1 ]->assemble_volume_integrals(cell);
        END_TIMER("assemble_volume_integrals");

        START_TIMER("assemble_fluxes_boundary");
        multidim_assembly_[ cell.dim()-1 ]->assemble_fluxes_boundary(cell);
        END_TIMER("assemble_fluxes_boundary");

        START_TIMER("assemble_fluxes_elem_elem");
        multidim_assembly_[ cell.dim()-1 ]->assemble_fluxes_element_element(cell);
        END_TIMER("assemble_fluxes_elem_elem");
    }

  START_TIMER("assemble_fluxes_elem_side");
    assemble_fluxes_element_side<1>();
    assemble_fluxes_element_side<2>();
    assemble_fluxes_element_side<3>();
  END_TIMER("assemble_fluxes_elem_side");
  END_TIMER("assemble_stiffness");
}



template<class Model>
void TransportDG<Model>::set_sources()
{
  START_TIMER("assemble_sources");
    Model::balance_->start_source_assembly(Model::subst_idx);
    set_sources<1>();
    set_sources<2>();
    set_sources<3>();
    Model::balance_->finish_source_assembly(Model::subst_idx);
  END_TIMER("assemble_sources");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_sources()
{
    FEValues<dim,3> fe_values(*assembly<dim>()->mapping(), *assembly<dim>()->quad(), *assembly<dim>()->fe(),
            update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = assembly<dim>()->fe()->n_dofs(), qsize = assembly<dim>()->quad()->size();
    vector<std::vector<double> > sources_conc(Model::n_substances(), std::vector<double>(qsize)),
            sources_density(Model::n_substances(), std::vector<double>(qsize)),
            sources_sigma(Model::n_substances(), std::vector<double>(qsize));
    vector<LongIdx> dof_indices(ndofs);
    vector<LongIdx> loc_dof_indices(ndofs);
    PetscScalar local_rhs[ndofs];
    vector<PetscScalar> local_source_balance_vector(ndofs), local_source_balance_rhs(ndofs);
    double source;

    // assemble integral over elements
    for (auto cell : data_->dh_->own_range() )
    {
        if (cell.dim() != dim) continue;
        ElementAccessor<3> elm = cell.elm();

        fe_values.reinit(elm);
        cell.get_dof_indices(dof_indices);
        cell.get_loc_dof_indices(loc_dof_indices);

        Model::compute_source_coefficients(fe_values.point_list(), elm, sources_conc, sources_density, sources_sigma);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
            fill_n(local_rhs, ndofs, 0);
            local_source_balance_vector.assign(ndofs, 0);
            local_source_balance_rhs.assign(ndofs, 0);

            // compute sources
            for (unsigned int k=0; k<qsize; k++)
            {
                source = (sources_density[sbi][k] + sources_conc[sbi][k]*sources_sigma[sbi][k])*fe_values.JxW(k);

                for (unsigned int i=0; i<ndofs; i++)
                    local_rhs[i] += source*fe_values.shape_value(i,k);
            }
            data_->ls[sbi]->rhs_set_values(ndofs, &(dof_indices[0]), local_rhs);

            for (unsigned int i=0; i<ndofs; i++)
            {
                for (unsigned int k=0; k<qsize; k++)
                    local_source_balance_vector[i] -= sources_sigma[sbi][k]*fe_values.shape_value(i,k)*fe_values.JxW(k);

                local_source_balance_rhs[i] += local_rhs[i];
            }
            Model::balance_->add_source_values(Model::subst_idx[sbi], elm.region().bulk_idx(), loc_dof_indices,
                                               local_source_balance_vector, local_source_balance_rhs);
        }
    }
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_element_side()
{

    if (dim == 1) return;
    FEValues<dim-1,3> fe_values_vb(*assembly<dim>()->mapping_low(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe_low(),
            update_values | update_gradients | update_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fe_values_side(*assembly<dim>()->mapping(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe(),
            update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*assembly<dim>()->mapping(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe_rt(),
            update_values);
    FEValues<dim-1,3> fv_rt(*assembly<dim>()->mapping_low(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe_rt_low(),
            update_values);

    vector<FEValuesSpaceBase<3>*> fv_sb(2);
    const unsigned int ndofs = assembly<dim>()->fe()->n_dofs();    // number of local dofs
    const unsigned int qsize = assembly<dim>()->quad_low()->size();     // number of quadrature points
    int side_dof_indices[2*ndofs];
    vector<LongIdx> indices(ndofs);
    unsigned int n_dofs[2], n_indices;
    vector<arma::vec3> velocity_higher, velocity_lower;
    vector<double> frac_sigma(qsize);
    vector<double> csection_lower(qsize), csection_higher(qsize);
    PetscScalar local_matrix[4*ndofs*ndofs];
    double comm_flux[2][2];

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    fv_sb[0] = &fe_values_vb;
    fv_sb[1] = &fe_values_side;

    // assemble integral over sides
    for (DHCellAccessor cell_lower_dim : data_->dh_->local_range() )
        for( DHCellSide neighb_side : cell_lower_dim.neighb_sides() )
        {
            // skip neighbours of different dimension
            if (cell_lower_dim.elm().dim() != dim-1) continue;

            ElementAccessor<3> elm_lower_dim = cell_lower_dim.elm();
            n_indices = cell_lower_dim.get_dof_indices(indices);
    		for(unsigned int i=0; i<n_indices; ++i) {
    			side_dof_indices[i] = indices[i];
    		}
            fe_values_vb.reinit(elm_lower_dim);
            n_dofs[0] = fv_sb[0]->n_dofs();

            DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.side()->element().idx() );
            ElementAccessor<3> elm_higher_dim = cell_higher_dim.elm();
            n_indices = cell_higher_dim.get_dof_indices(indices);
    		for(unsigned int i=0; i<n_indices; ++i) {
    			side_dof_indices[i+n_dofs[0]] = indices[i];
    		}
            fe_values_side.reinit(elm_higher_dim, neighb_side.side()->side_idx());
            n_dofs[1] = fv_sb[1]->n_dofs();

            // Testing element if they belong to local partition.
            bool own_element_id[2];
            own_element_id[0] = cell_lower_dim.is_own();
            own_element_id[1] = cell_higher_dim.is_own();

            fsv_rt.reinit(elm_higher_dim, neighb_side.side()->side_idx());
            fv_rt.reinit(elm_lower_dim);
            calculate_velocity(elm_higher_dim, velocity_higher, fsv_rt);
            calculate_velocity(elm_lower_dim, velocity_lower, fv_rt);
            Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_lower, elm_lower_dim, data_->ad_coef_edg[0], data_->dif_coef_edg[0]);
            Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_higher, elm_higher_dim, data_->ad_coef_edg[1], data_->dif_coef_edg[1]);
            data_->cross_section.value_list(fe_values_vb.point_list(), elm_lower_dim, csection_lower);
            data_->cross_section.value_list(fe_values_vb.point_list(), elm_higher_dim, csection_higher);

            for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++) // Optimize: SWAP LOOPS
            {
                for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
                    for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
                        local_matrix[i*(n_dofs[0]+n_dofs[1])+j] = 0;

                data_->fracture_sigma[sbi].value_list(fe_values_vb.point_list(), elm_lower_dim, frac_sigma);

                // set transmission conditions
                for (unsigned int k=0; k<qsize; k++)
                {
                    /* The communication flux has two parts:
                    * - "diffusive" term containing sigma
                    * - "advective" term representing usual upwind
                    *
                    * The calculation differs from the reference manual, since ad_coef and dif_coef have different meaning
                    * than b and A in the manual.
                    * In calculation of sigma there appears one more csection_lower in the denominator.
                    */
                    double sigma = frac_sigma[k]*arma::dot(data_->dif_coef_edg[0][sbi][k]*fe_values_side.normal_vector(k),fe_values_side.normal_vector(k))*
                            2*csection_higher[k]*csection_higher[k]/(csection_lower[k]*csection_lower[k]);

                    double transport_flux = arma::dot(data_->ad_coef_edg[1][sbi][k], fe_values_side.normal_vector(k));

                    comm_flux[0][0] =  (sigma-min(0.,transport_flux))*fv_sb[0]->JxW(k);
                    comm_flux[0][1] = -(sigma-min(0.,transport_flux))*fv_sb[0]->JxW(k);
                    comm_flux[1][0] = -(sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);
                    comm_flux[1][1] =  (sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);

                    for (int n=0; n<2; n++)
                    {
                        if (!own_element_id[n]) continue;

                        for (unsigned int i=0; i<n_dofs[n]; i++)
                            for (int m=0; m<2; m++)
                                for (unsigned int j=0; j<n_dofs[m]; j++)
                                    local_matrix[(i+n*n_dofs[0])*(n_dofs[0]+n_dofs[1]) + m*n_dofs[0] + j] +=
                                            comm_flux[m][n]*fv_sb[m]->shape_value(j,k)*fv_sb[n]->shape_value(i,k);
                    }
                }
                data_->ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], side_dof_indices, n_dofs[0]+n_dofs[1], side_dof_indices, local_matrix);
            }
        }

}






template<class Model>
void TransportDG<Model>::set_boundary_conditions()
{
  START_TIMER("assemble_bc");
    Model::balance_->start_flux_assembly(Model::subst_idx);
    set_boundary_conditions<1>();
    set_boundary_conditions<2>();
    set_boundary_conditions<3>();
    Model::balance_->finish_flux_assembly(Model::subst_idx);
  END_TIMER("assemble_bc");
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_boundary_conditions()
{
    FESideValues<dim,3> fe_values_side(*assembly<dim>()->mapping(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe(),
            update_values | update_gradients | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*assembly<dim>()->mapping(), *assembly<dim>()->quad_low(), *assembly<dim>()->fe_rt(),
                update_values);
    const unsigned int ndofs = assembly<dim>()->fe()->n_dofs(), qsize = assembly<dim>()->quad_low()->size();
    vector<LongIdx> side_dof_indices(ndofs);
    unsigned int loc_b=0;
    double local_rhs[ndofs];
    vector<PetscScalar> local_flux_balance_vector(ndofs);
    PetscScalar local_flux_balance_rhs;
    vector<double> bc_values(qsize);
    vector<double> bc_fluxes(qsize),
            bc_sigma(qsize),
            bc_ref_values(qsize);
    vector<double> csection(qsize);
    vector<arma::vec3> velocity;

    for (auto cell : data_->dh_->own_range() )
    {
        if (cell.elm()->boundary_idx_ == nullptr) continue;

        for (unsigned int si=0; si<cell.elm()->n_sides(); si++)
        {
            const Edge *edg = cell.elm().side(si)->edge();
            if (edg->n_sides > 1) continue;
            // skip edges lying not on the boundary
            if (edg->side(0)->cond() == NULL) continue;

            if (edg->side(0)->dim() != dim-1)
            {
                if (edg->side(0)->cond() != nullptr) ++loc_b;
                continue;
            }

            SideIter side = edg->side(0);
            ElementAccessor<3> elm = Model::mesh_->element_accessor( side->element().idx() );
            ElementAccessor<3> ele_acc = side->cond()->element_accessor();

            arma::uvec bc_type;
            Model::get_bc_type(ele_acc, bc_type);

            fe_values_side.reinit(elm, side->side_idx());
            fsv_rt.reinit(elm, side->side_idx());
            calculate_velocity(elm, velocity, fsv_rt);

            Model::compute_advection_diffusion_coefficients(fe_values_side.point_list(), velocity, side->element(), data_->ad_coef, data_->dif_coef);
            data_->cross_section.value_list(fe_values_side.point_list(), side->element(), csection);

            auto dh_cell = data_->dh_->cell_accessor_from_element( side->element().idx() );
            dh_cell.get_dof_indices(side_dof_indices);

            for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
            {
                fill_n(local_rhs, ndofs, 0);
                local_flux_balance_vector.assign(ndofs, 0);
                local_flux_balance_rhs = 0;
                
                // The b.c. data are fetched for all possible b.c. types since we allow
                // different bc_type for each substance.
                data_->bc_dirichlet_value[sbi].value_list(fe_values_side.point_list(), ele_acc, bc_values);

                double side_flux = 0;
                for (unsigned int k=0; k<qsize; k++)
                    side_flux += arma::dot(data_->ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
                double transport_flux = side_flux/side->measure();

                if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux < 0)
                {
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        double bc_term = -transport_flux*bc_values[k]*fe_values_side.JxW(k);
                        for (unsigned int i=0; i<ndofs; i++)
                            local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
                    }
                    for (unsigned int i=0; i<ndofs; i++)
                        local_flux_balance_rhs -= local_rhs[i];
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
                {
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        double bc_term = data_->gamma[sbi][side->cond_idx()]*bc_values[k]*fe_values_side.JxW(k);
                        arma::vec3 bc_grad = -bc_values[k]*fe_values_side.JxW(k)*data_->dg_variant*(arma::trans(data_->dif_coef[sbi][k])*fe_values_side.normal_vector(k));
                        for (unsigned int i=0; i<ndofs; i++)
                            local_rhs[i] += bc_term*fe_values_side.shape_value(i,k)
                                    + arma::dot(bc_grad,fe_values_side.shape_grad(i,k));
                    }
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        for (unsigned int i=0; i<ndofs; i++)
                        {
                            local_flux_balance_vector[i] += (arma::dot(data_->ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.shape_value(i,k)
                                    - arma::dot(data_->dif_coef[sbi][k]*fe_values_side.shape_grad(i,k),fe_values_side.normal_vector(k))
                                    + data_->gamma[sbi][side->cond_idx()]*fe_values_side.shape_value(i,k))*fe_values_side.JxW(k);
                        }
                    }
                    if (Model::time_->tlevel() > 0)
                        for (unsigned int i=0; i<ndofs; i++)
                            local_flux_balance_rhs -= local_rhs[i];
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
                {
                    Model::get_flux_bc_data(sbi, fe_values_side.point_list(), ele_acc, bc_fluxes, bc_sigma, bc_ref_values);
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        double bc_term = csection[k]*(bc_sigma[k]*bc_ref_values[k]+bc_fluxes[k])*fe_values_side.JxW(k);
                        for (unsigned int i=0; i<ndofs; i++)
                            local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
                    }

                    for (unsigned int i=0; i<ndofs; i++)
                    {
                        for (unsigned int k=0; k<qsize; k++)
                            local_flux_balance_vector[i] += csection[k]*bc_sigma[k]*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
                        local_flux_balance_rhs -= local_rhs[i];
                    }
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
                {
                    Model::get_flux_bc_data(sbi, fe_values_side.point_list(), ele_acc, bc_fluxes, bc_sigma, bc_ref_values);
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        double bc_term = csection[k]*(bc_sigma[k]*bc_ref_values[k]+bc_fluxes[k])*fe_values_side.JxW(k);
                        for (unsigned int i=0; i<ndofs; i++)
                            local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
                    }

                    for (unsigned int i=0; i<ndofs; i++)
                    {
                        for (unsigned int k=0; k<qsize; k++)
                            local_flux_balance_vector[i] += csection[k]*(arma::dot(data_->ad_coef[sbi][k], fe_values_side.normal_vector(k)) + bc_sigma[k])*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
                        local_flux_balance_rhs -= local_rhs[i];
                    }
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux >= 0)
                {
                    for (unsigned int k=0; k<qsize; k++)
                    {
                        for (unsigned int i=0; i<ndofs; i++)
                            local_flux_balance_vector[i] += arma::dot(data_->ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
                    }
                }
                data_->ls[sbi]->rhs_set_values(ndofs, &(side_dof_indices[0]), local_rhs);

                Model::balance_->add_flux_matrix_values(Model::subst_idx[sbi], loc_b, side_dof_indices, local_flux_balance_vector);
                Model::balance_->add_flux_vec_value(Model::subst_idx[sbi], loc_b, local_flux_balance_rhs);
            }
            ++loc_b;
        }
    }
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::calculate_velocity(const ElementAccessor<3> &cell, 
                                            vector<arma::vec3> &velocity, 
                                            FEValuesBase<dim,3> &fv)
{
    ASSERT_EQ(cell->dim(), dim).error("Element dimension mismatch!");

    velocity.resize(fv.n_points());
    arma::mat map_mat = assembly<dim>()->mapping()->element_map(cell);
    vector<arma::vec3> point_list;
    point_list.resize(fv.n_points());
    for (unsigned int k=0; k<fv.n_points(); k++)
    	point_list[k] = assembly<dim>()->mapping()->project_unit_to_real(RefElement<dim>::local_to_bary(fv.get_quadrature()->point(k)), map_mat);
    Model::velocity_field_ptr_->value_list(point_list, cell, velocity);
}



template<class Model>
void TransportDG<Model>::set_initial_condition()
{
    START_TIMER("set_init_cond");
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        data_->ls[sbi]->start_allocation();
    prepare_initial_condition<1>();
    prepare_initial_condition<2>();
    prepare_initial_condition<3>();

    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        data_->ls[sbi]->start_add_assembly();
    prepare_initial_condition<1>();
    prepare_initial_condition<2>();
    prepare_initial_condition<3>();

    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    {
        data_->ls[sbi]->finish_assembly();
        data_->ls[sbi]->solve();
    }
    END_TIMER("set_init_cond");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::prepare_initial_condition()
{
    FEValues<dim,3> fe_values(*assembly<dim>()->mapping(), *assembly<dim>()->quad(), *assembly<dim>()->fe(),
            update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = assembly<dim>()->fe()->n_dofs(), qsize = assembly<dim>()->quad()->size();
    std::vector<LongIdx> dof_indices(ndofs);
    double matrix[ndofs*ndofs], rhs[ndofs];
    std::vector<std::vector<double> > init_values(Model::n_substances());

    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++) // Optimize: SWAP LOOPS
        init_values[sbi].resize(qsize);

    for (auto cell : data_->dh_->own_range() )
    {
        if (cell.dim() != dim) continue;
        ElementAccessor<3> elem = cell.elm();

        cell.get_dof_indices(dof_indices);
        fe_values.reinit(elem);

        Model::compute_init_cond(fe_values.point_list(), elem, init_values);

        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
            for (unsigned int i=0; i<ndofs; i++)
            {
                rhs[i] = 0;
                for (unsigned int j=0; j<ndofs; j++)
                    matrix[i*ndofs+j] = 0;
            }

            for (unsigned int k=0; k<qsize; k++)
            {
                double rhs_term = init_values[sbi][k]*fe_values.JxW(k);

                for (unsigned int i=0; i<ndofs; i++)
                {
                    for (unsigned int j=0; j<ndofs; j++)
                        matrix[i*ndofs+j] += fe_values.shape_value(i,k)*fe_values.shape_value(j,k)*fe_values.JxW(k);

                    rhs[i] += fe_values.shape_value(i,k)*rhs_term;
                }
            }
            data_->ls[sbi]->set_values(ndofs, &(dof_indices[0]), ndofs, &(dof_indices[0]), matrix, rhs);
        }
    }
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

            unsigned int n_dofs;
            switch (cell.dim())
            {
            case 1:
                n_dofs = assembly1_->fe()->n_dofs();
                break;
            case 2:
                n_dofs = assembly2_->fe()->n_dofs();
                break;
            case 3:
                n_dofs = assembly3_->fe()->n_dofs();
                break;
            }

			std::vector<LongIdx> dof_indices(n_dofs);
            cell.get_dof_indices(dof_indices);

            for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
            {
                double old_average = 0;
                for (unsigned int j=0; j<n_dofs; ++j)
                    old_average += data_->ls[sbi]->get_solution_array()[dof_indices[j]-data_->dh_->distr()->begin()];
                old_average /= n_dofs;

                for (unsigned int j=0; j<n_dofs; ++j)
                    data_->ls[sbi]->get_solution_array()[dof_indices[j]-data_->dh_->distr()->begin()] += solution_elem_[sbi][i_cell] - old_average;
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

template<class Model>
void TransportDG<Model>::initialize_assembly_objects()
{
    for (unsigned int i=0; i<multidim_assembly_.size(); ++i)
    	multidim_assembly_[i]->initialize();
}

template<class Model>
template<unsigned int dim>
std::shared_ptr<AssemblyDG<dim, Model>> TransportDG<Model>::assembly() {
    return std::dynamic_pointer_cast<AssemblyDG<dim, Model>>(multidim_assembly_[dim-1]);
}







template class TransportDG<ConcentrationTransportModel>;
template class TransportDG<HeatTransferModel>;




