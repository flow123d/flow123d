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
#include "mechanics/elasticity.hh"
#include "mechanics/assembly_elasticity.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_system.hh"
#include "fields/field_fe.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_PERMON.hh"
#include "coupling/balance.hh"
#include "mesh/neighbours.h"
#include "coupling/generic_assembly.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "fields/field_model.hh"
#include "input/factory.hh"




using namespace Input::Type;



const Record & Elasticity::get_input_type() {
    std::string equation_name = std::string(name_) + "_FE";
	return IT::Record(
                std::string(equation_name),
                "FEM for linear elasticity.")
           .copy_keys(EquationBase::record_template())
		   .copy_keys(EquationBase::user_fields_template(equation_name))
           .declare_key("balance", Balance::get_input_type(), Default("{}"),
                    "Settings for computing balance.")
           .declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
                    "Parameters of output stream.")
           .declare_key("solver", LinSys_PETSC::get_input_type(), Default::obligatory(),
				"Linear solver for elasticity.")
		   .declare_key("input_fields", Array(
		        Elasticity::EqFields()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the equation.")
           .declare_key("output",
                EqFields().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ \"displacement\" ] }"),
                "Setting of the field output.")
           .declare_key("contact", Bool(), IT::Default("false"), "Indicates the use of contact conditions on fractures.")
		   .close();
}

const int Elasticity::registrar =
		Input::register_class< Elasticity, Mesh &, const Input::Record>(std::string(name_) + "_FE") +
		Elasticity::get_input_type().size();



double lame_mu(double young, double poisson)
{
    return young*0.5/(poisson+1.);
}


double lame_lambda(double young, double poisson)
{
    return young*poisson/((poisson+1.)*(1.-2.*poisson));
}

// Functor computing lame_mu
struct fn_lame_mu {
	inline double operator() (double young, double poisson) {
        return young * 0.5 / (poisson+1.);
    }
};

// Functor computing lame_lambda
struct fn_lame_lambda {
	inline double operator() (double young, double poisson) {
        return young * poisson / ((poisson+1.)*(1.-2.*poisson));
    }
};

// Functor computing base of dirichlet_penalty (without dividing by side meassure)
struct fn_dirichlet_penalty {
	inline double operator() (double lame_mu, double lame_lambda) {
        return 1e3 * (2 * lame_mu + lame_lambda);
    }
};






const Selection & Elasticity::EqFields::get_bc_type_selection() {
    return Selection("Elasticity_BC_Type", "Types of boundary conditions for mechanics.")
            .add_value(bc_type_displacement, "displacement",
                  "Prescribed displacement.")
            .add_value(bc_type_displacement_normal, "displacement_n",
                  "Prescribed displacement in the normal direction to the boundary.")
            .add_value(bc_type_traction, "traction",
                  "Prescribed traction.")
            .add_value(bc_type_stress, "stress",
                  "Prescribed stress tensor.")
            .close();
}


Elasticity::EqFields::EqFields()
{
    *this+=bc_type
        .name("bc_type")
        .description(
        "Type of boundary condition.")
        .units( UnitSI::dimensionless() )
        .input_default("\"traction\"")
        .input_selection( get_bc_type_selection() )
        .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
        
    *this+=bc_displacement
        .name("bc_displacement")
        .description("Prescribed displacement on boundary.")
        .units( UnitSI().m() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=bc_traction
        .name("bc_traction")
        .description("Prescribed traction on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=bc_stress
        .name("bc_stress")
        .description("Prescribed stress on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=load
        .name("load")
        .description("Prescribed bulk load.")
        .units( UnitSI().N().m(-3) )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=young_modulus
        .name("young_modulus")
        .description("Young's modulus.")
        .units( UnitSI().Pa() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
    
    *this+=poisson_ratio
        .name("poisson_ratio")
        .description("Poisson's ratio.")
        .units( UnitSI().dimensionless() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
        
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of transfer of forces through fractures.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_rhs);

    *this+=initial_stress
        .name("initial_stress")
        .description("Initial stress tensor.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);

    *this += region_id.name("region_id")
    	        .units( UnitSI::dimensionless())
    	        .flags(FieldFlag::equation_external_output);
                
    *this += subdomain.name("subdomain")
      .units( UnitSI::dimensionless() )
      .flags(FieldFlag::equation_external_output);
      
    *this+=cross_section
      .name("cross_section")
      .units( UnitSI().m(3).md() )
      .flags(input_copy & in_time_term & in_main_matrix & in_rhs);
    
    *this+=cross_section_min
      .name("cross_section_min")
      .description("Minimal cross-section of fractures.")
      .units( UnitSI().m(3).md() )
      .input_default("0.0");
      
    *this+=potential_load
      .name("potential_load")
      .units( UnitSI().m() )
      .flags(input_copy & in_rhs);

    *this+=output_field
            .name("displacement")
            .description("Displacement vector field output.")
            .units( UnitSI().m() )
            .flags(equation_result);
    
    *this += output_stress
            .name("stress")
            .description("Stress tensor output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_von_mises_stress
            .name("von_mises_stress")
            .description("von Mises stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_mean_stress
            .name("mean_stress")
            .description("mean stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);

    *this += output_cross_section
            .name("cross_section_updated")
            .description("Cross-section after deformation - output.")
            .units( UnitSI().m() )
            .flags(equation_result);
            
    *this += output_divergence
            .name("displacement_divergence")
            .description("Displacement divergence output.")
            .units( UnitSI().dimensionless() )
            .flags(equation_result);

    *this += lame_mu.name("lame_mu")
            .description("Field lame_mu.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += lame_lambda.name("lame_lambda")
            .description("Field lame_lambda.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += dirichlet_penalty.name("dirichlet_penalty")
            .description("Field dirichlet_penalty.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    // add all input fields to the output list
    output_fields += *this;

    this->add_coords_field();
    this->set_default_fieldset();

}

void Elasticity::EqData::create_dh(Mesh * mesh, unsigned int fe_order)
{
	ASSERT_EQ(fe_order, 1)(fe_order).error("Unsupported polynomial order for finite elements in Elasticity");
    MixedPtr<FE_P> fe_p(1);
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);

    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
	dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh);

	dh_->distribute_dofs(ds);


    MixedPtr<FE_P_disc> fe_p_disc(0);
    dh_scalar_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> ds_scalar = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_p_disc);
	dh_scalar_->distribute_dofs(ds_scalar);


    MixedPtr<FiniteElement> fe_t = mixed_fe_system(MixedPtr<FE_P_disc>(0), FEType::FETensor, 9);
    dh_tensor_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> dst = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_t);
	dh_tensor_->distribute_dofs(dst);
}


Elasticity::Elasticity(Mesh & init_mesh, const Input::Record in_rec, TimeGovernor *tm)
        : EquationBase(init_mesh, in_rec),
		  input_rec(in_rec),
		  stiffness_assembly_(nullptr),
		  rhs_assembly_(nullptr),
          constraint_assembly_(nullptr),
          cs_assembly_(nullptr),
		  output_fields_assembly_(nullptr)
{
	// Can not use name() + "constructor" here, since START_TIMER only accepts const char *
	// due to constexpr optimization.
	START_TIMER(name_);

    eq_data_ = std::make_shared<EqData>();
    eq_fields_ = std::make_shared<EqFields>();
    this->eq_fieldset_ = eq_fields_;
    
    auto time_rec = in_rec.val<Input::Record>("time");
    if (tm == nullptr)
    {
        time_ = new TimeGovernor(time_rec);
    }
    else
    {
        TimeGovernor time_from_rec(time_rec);
        ASSERT( time_from_rec.is_default() ).error("Duplicate key 'time', time in elasticity is already initialized from parent class!");
        time_ = tm;
    }


    // Set up physical parameters.
    eq_fields_->set_mesh(init_mesh);
    eq_fields_->region_id = GenericField<3>::region_id(*mesh_);
    eq_fields_->subdomain = GenericField<3>::subdomain(*mesh_);
    eq_data_->balance_ = this->balance();
    
    // create finite element structures and distribute DOFs
    eq_data_->create_dh(mesh_, 1);
    DebugOut().fmt("Mechanics: solution size {}\n", eq_data_->dh_->n_global_dofs());
    
}


void Elasticity::initialize()
{
    output_stream_ = OutputTime::create_output_stream("mechanics", input_rec.val<Input::Record>("output_stream"), time().get_unit_conversion());
    
    eq_fields_->set_components({"displacement"});
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), time() );
    
//     balance_ = std::make_shared<Balance>("mechanics", mesh_);
//     balance_->init_from_input(input_rec.val<Input::Record>("balance"), *time_);
    // initialization of balance object
//     eq_data_->balance_idx_ = {
//         balance_->add_quantity("force_x"),
//         balance_->add_quantity("force_y"),
//         balance_->add_quantity("force_z")
//     };
//     balance_->units(UnitSI().kg().m().s(-2));

    // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
    eq_fields_->output_field_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_);
    eq_fields_->output_field.set(eq_fields_->output_field_ptr, 0.);
    
    // setup output stress
    eq_fields_->output_stress_ptr = create_field_fe<3, FieldValue<3>::TensorFixed>(eq_data_->dh_tensor_);
    eq_fields_->output_stress.set(eq_fields_->output_stress_ptr, 0.);
    
    // setup output von Mises stress
    eq_fields_->output_von_mises_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_von_mises_stress.set(eq_fields_->output_von_mises_stress_ptr, 0.);

    // setup output mean stress
    eq_fields_->output_mean_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_mean_stress.set(eq_fields_->output_mean_stress_ptr, 0.);
    
    // setup output cross-section
    eq_fields_->output_cross_section_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_cross_section.set(eq_fields_->output_cross_section_ptr, 0.);
    
    // setup output divergence
    eq_fields_->output_div_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_divergence.set(eq_fields_->output_div_ptr, 0.);

    // read optional user fields
    Input::Array user_fields_arr;
    if (input_rec.opt_val("user_fields", user_fields_arr)) {
       	this->init_user_fields(user_fields_arr, eq_fields_->output_fields);
    }
    
    eq_fields_->output_fields.set_mesh(*mesh_);
    eq_fields_->output_field.output_type(OutputTime::CORNER_DATA);

    // set time marks for writing the output
    eq_fields_->output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), this->time());

    // set instances of FieldModel
    eq_fields_->lame_mu.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_mu(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->lame_lambda.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_lambda(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->dirichlet_penalty.set(Model<3, FieldValue<3>::Scalar>::create(fn_dirichlet_penalty(), eq_fields_->lame_mu, eq_fields_->lame_lambda), 0.0);

    // equation default PETSc solver options
    std::string petsc_default_opts;
    petsc_default_opts = "-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg";
    
    // allocate matrix and vector structures
    LinSys *ls;
    has_contact_ = input_rec.val<bool>("contact");
    if (has_contact_) {
#ifndef FLOW123D_HAVE_PERMON
        ASSERT(false).error("Flow123d was not built with PERMON library, therefore contact conditions are unsupported.");
#endif //FLOW123D_HAVE_PERMON
        ls = new LinSys_PERMON(eq_data_->dh_->distr().get(), petsc_default_opts);

        // allocate constraint matrix and vector
        unsigned int n_own_constraints = 0; // count locally owned cells with neighbours
        for (auto cell : eq_data_->dh_->own_range())
            if (cell.elm()->n_neighs_vb() > 0)
                n_own_constraints++;
        unsigned int n_constraints = 0; // count all cells with neighbours
        for (auto elm : mesh_->elements_range())
            if (elm->n_neighs_vb() > 0)
                eq_data_->constraint_idx[elm.idx()] = n_constraints++;
        unsigned int nnz = eq_data_->dh_->ds()->fe()[1_d]->n_dofs()*mesh_->max_edge_sides(1) +
                        eq_data_->dh_->ds()->fe()[2_d]->n_dofs()*mesh_->max_edge_sides(2) +
                        eq_data_->dh_->ds()->fe()[3_d]->n_dofs()*mesh_->max_edge_sides(3);
        MatCreateAIJ(PETSC_COMM_WORLD, n_own_constraints, eq_data_->dh_->lsize(), PETSC_DECIDE, PETSC_DECIDE, nnz, 0, nnz, 0, &eq_data_->constraint_matrix);
        VecCreateMPI(PETSC_COMM_WORLD, n_own_constraints, PETSC_DECIDE, &eq_data_->constraint_vec);
        ((LinSys_PERMON*)ls)->set_inequality(eq_data_->constraint_matrix,eq_data_->constraint_vec);

        constraint_assembly_ = new GenericAssembly< ConstraintAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    } else {
        ls = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ((LinSys_PETSC*)ls)->set_initial_guess_nonzero();
    }
    ls->set_from_input( input_rec.val<Input::Record>("solver") );
    ls->set_solution(eq_fields_->output_field_ptr->vec().petsc_vec());
    eq_data_->ls = ls;

    stiffness_assembly_ = new GenericAssembly< StiffnessAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    rhs_assembly_ = new GenericAssembly< RhsAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    cs_assembly_ = new GenericAssembly< CrossSectionAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    output_fields_assembly_ = new GenericAssembly< OutpuFieldsAssemblyElasticity >(eq_fields_.get(), eq_data_.get());

    // initialization of balance object
//     balance_->allocate(eq_data_->dh_->distr()->lsize(),
//             max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs())));

}


Elasticity::~Elasticity()
{
//     delete time_;

    if (stiffness_assembly_!=nullptr) delete stiffness_assembly_;
    if (rhs_assembly_!=nullptr) delete rhs_assembly_;
    if (constraint_assembly_ != nullptr) delete constraint_assembly_;
    if (cs_assembly_!=nullptr) delete cs_assembly_;
    if (output_fields_assembly_!=nullptr) delete output_fields_assembly_;

    eq_data_.reset();
    eq_fields_.reset();
}



void Elasticity::update_output_fields()
{
    eq_fields_->set_time(time_->step(), LimitSide::right);
    
    // update ghost values of solution vector and prepare dependent fields
	eq_fields_->output_field_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().zero_entries();
	eq_fields_->output_cross_section_ptr->vec().zero_entries();
	eq_fields_->output_div_ptr->vec().zero_entries();
	eq_fields_->output_field_ptr->vec().local_to_ghost_end();

    
    // compute first the updated cross-section
    cs_assembly_->assemble(eq_data_->dh_);
    // then compute output fields depending on solution and/or on cross-section (stress, divergence etc.)
    output_fields_assembly_->assemble(eq_data_->dh_);

    // update ghost values of computed fields
    eq_fields_->output_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_div_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_end();
    eq_fields_->output_div_ptr->vec().local_to_ghost_end();
}




void Elasticity::zero_time_step()
{
	START_TIMER(name_);
	eq_fields_->mark_input_times( *time_ );
	eq_fields_->set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "mechanics") ) {
		WarningOut() << ss.str();
	}

    preallocate();
    

    // after preallocation we assemble the matrices and vectors required for balance of forces
//     for (auto subst_idx : eq_data_->balance_idx_)
//         balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
    
//     update_solution();
    eq_data_->ls->start_add_assembly();
    MatSetOption(*eq_data_->ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    eq_data_->ls->mat_zero_entries();
    eq_data_->ls->rhs_zero_entries();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);
    eq_data_->ls->finish_assembly();
    LinSys::SolveInfo si = eq_data_->ls->solve();
    MessageOut().fmt("[mech solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, eq_data_->ls->compute_residual());
    output_data();
}



void Elasticity::preallocate()
{
    // preallocate system matrix
	eq_data_->ls->start_allocation();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);

    if (has_contact_)
        assemble_constraint_matrix();
}




void Elasticity::update_solution()
{
	START_TIMER("DG-ONE STEP");

    next_time();
	solve_linear_system();

    calculate_cumulative_balance();
    
    output_data();

    END_TIMER("DG-ONE STEP");
}

void Elasticity::next_time()
{
    time_->next_time();
    time_->view("MECH");
    
}



void Elasticity::solve_linear_system()
{
    START_TIMER("data reinit");
    eq_fields_->set_time(time_->step(), LimitSide::right);
    END_TIMER("data reinit");
    
    // assemble stiffness matrix
    if (eq_data_->ls->get_matrix() == NULL
        || eq_fields_->subset(FieldFlag::in_main_matrix).changed())
    {
        DebugOut() << "Mechanics: Assembling matrix.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->mat_zero_entries();
        stiffness_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (eq_data_->ls->get_rhs() == NULL
        || eq_fields_->subset(FieldFlag::in_rhs).changed())
    {
        DebugOut() << "Mechanics: Assembling right hand side.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->rhs_zero_entries();
        rhs_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

    START_TIMER("solve");
    LinSys::SolveInfo si = eq_data_->ls->solve();
    MessageOut().fmt("[mech solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, eq_data_->ls->compute_residual());
    END_TIMER("solve");
}






void Elasticity::output_data()
{
    START_TIMER("MECH-OUTPUT");

    // gather the solution from all processors
    eq_fields_->output_fields.set_time( this->time().step(), LimitSide::left);
    //if (eq_fields_->output_fields.is_field_output_time(eq_fields_->output_field, this->time().step()) )
    update_output_fields();
    eq_fields_->output_fields.output(this->time().step());

//     START_TIMER("MECH-balance");
//     balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
//     balance_->output();
//     END_TIMER("MECH-balance");

    END_TIMER("MECH-OUTPUT");
}




void Elasticity::calculate_cumulative_balance()
{
//     if (balance_->cumulative())
//     {
//         balance_->calculate_cumulative(subst_idx, eq_data_->ls->get_solution());
//     }
}


void Elasticity::assemble_constraint_matrix()
{
    MatZeroEntries(eq_data_->constraint_matrix);
    VecZeroEntries(eq_data_->constraint_vec);
    constraint_assembly_->assemble(eq_data_->dh_);
    MatAssemblyBegin(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(eq_data_->constraint_vec);
    VecAssemblyEnd(eq_data_->constraint_vec);
}


