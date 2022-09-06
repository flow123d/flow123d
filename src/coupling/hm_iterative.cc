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
 * @file    hm_iterative.cc
 * @brief   
 * @author  Jan Stebel
 */

#include "hm_iterative.hh"
#include "system/sys_profiler.hh"
#include "input/input_type.hh"
#include "flow/darcy_flow_lmh.hh"
#include "fields/field_fe.hh"         // for create_field_fe()
#include "fields/field_model.hh"      // for Model
#include "assembly_hm.hh"


FLOW123D_FORCE_LINK_IN_CHILD(coupling_iterative)


namespace it = Input::Type;


struct fn_pressure_potential {
    inline double operator() (double alpha, double density, double gravity, double pressure)
    {
        return -alpha*density*gravity*pressure;
    }
};


struct fn_hm_coupling_beta {

    fn_hm_coupling_beta(double beta_f) : beta_factor(beta_f) {}


    inline double operator() (double alpha, double lame_mu, double lame_lambda, double density, double gravity)
    {
        return beta_factor*0.5*alpha*alpha/(2*lame_mu/3 + lame_lambda)*density*gravity;
    }

private:

    const double beta_factor; ///< User-defined factor for iteration parameter.

};


struct fn_fluid_source {

    fn_fluid_source(const TimeGovernor *time) : time_(time) {}

    inline double operator() (double alpha, double beta, double pressure, double old_pressure, double div_u, double old_div_u)
    {
        return (beta*(pressure-old_pressure) + alpha*(old_div_u - div_u)) / time_->dt();
    }

private:

    const TimeGovernor *time_;
};


const it::Record & HM_Iterative::get_input_type() {
    std::string equation_name = "Coupling_Iterative";
    return it::Record(equation_name,
            "Record with data for iterative coupling of flow and mechanics.\n")
        .derive_from( DarcyFlowInterface::get_input_type() )
        .copy_keys(EquationBase::record_template())
		.copy_keys(EquationBase::user_fields_template(equation_name))
        .copy_keys(IterativeCoupling::record_template())
		.declare_key("flow_equation", DarcyLMH::get_input_type(),
		        it::Default::obligatory(),
				"Flow equation, provides the velocity field as a result.")
		.declare_key("mechanics_equation", Elasticity::get_input_type(),
				it::Default::obligatory(),
				"Mechanics, provides the displacement field.")
        .declare_key("input_fields", it::Array(
		        HM_Iterative::EqFields()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the HM coupling.")
        .declare_key( "iteration_parameter", it::Double(), it::Default("1"),
                "Tuning parameter for iterative splitting. Its default value "
                "corresponds to a theoretically optimal value with fastest convergence." )
        .declare_key( "max_it", it::Integer(0), it::Default("100"),
                "Maximal count of HM iterations." )
        .declare_key( "min_it", it::Integer(0), it::Default("1"),
                "Minimal count of HM iterations." )
        .declare_key( "a_tol", it::Double(0), it::Default("0"),
                "Absolute tolerance for difference in HM iteration." )
        .declare_key( "r_tol", it::Double(0), it::Default("1e-7"),
                "Relative tolerance for difference in HM iteration." )
		.close();
}


const int HM_Iterative::registrar = Input::register_class< HM_Iterative, Mesh &, const Input::Record >("Coupling_Iterative")
                                    + HM_Iterative::get_input_type().size();


HM_Iterative::EqFields::EqFields()
{
    *this += alpha.name("biot_alpha")
                     .description("Biot poroelastic coefficient.")
                     .units(UnitSI().dimensionless())
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += density.name("fluid_density")
                     .description("Volumetric mass density of the fluid.")
                     .units(UnitSI().kg().m(-3))
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += gravity.name("gravity")
                     .description("Gravitational acceleration constant.")
                     .units(UnitSI().m().s(-2))
                     .input_default("9.81")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += beta.name("relaxation_beta")
                     .description("Parameter of numerical method for iterative solution of hydro-mechanical coupling.")
                     .units(UnitSI().dimensionless())
                     .flags(FieldFlag::equation_external_output);
    
    *this += pressure_potential.name("pressure_potential")
                     .description("Coupling term entering the mechanics equation.")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_result);
    
    *this += old_pressure.name("old_pressure")
                     .description("Pressure from last computed time.")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_external_output);
    
    *this += old_iter_pressure.name("old_iter_pressure")
                     .description("Pressure from last computed iteration.")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_external_output);

    *this += old_div_u.name("old_displacement_divergence")
                     .description("Displacement divergence from last computed time.")
                     .units(UnitSI().dimensionless())
                     .flags(FieldFlag::equation_external_output);

    *this += ref_pressure_potential.name("ref_pressure_potential")
                     .description("Pressure potential on boundary (taking into account the flow boundary condition.")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_result);
    
    *this += flow_source.name("extra_flow_source")
                     .description("Coupling term entering the flow equation.")
                     .units(UnitSI().s(-1))
                     .flags(FieldFlag::equation_result);

    this->set_default_fieldset();
}


void HM_Iterative::EqFields::initialize(Mesh &mesh, HM_Iterative::EqData &eq_data, const TimeGovernor *time_, double beta_)
{
    // initialize coupling fields with FieldFE
    set_mesh(mesh);

    pressure_potential.set(Model<3, FieldValue<3>::Scalar>::create(
        fn_pressure_potential(),
        alpha,
        density,
        gravity,
        eq_data.flow_->eq_fields().field_edge_pressure
        ), 0.0);
    
    ref_potential_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(mesh, MixedPtr<FE_CR>());
    ref_pressure_potential.set(ref_potential_ptr_, 0.0);

    beta.set(Model<3, FieldValue<3>::Scalar>::create(
        fn_hm_coupling_beta(beta_),
        alpha,
        eq_data.mechanics_->eq_fields().lame_mu,
        eq_data.mechanics_->eq_fields().lame_lambda,
        density,
        gravity
        ), 0.0);

    auto old_pressure_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(eq_data.flow_->eq_data().dh_cr_, &eq_data.flow_->eq_data().p_edge_solution_previous_time);
    old_pressure.set(old_pressure_ptr_, 0.0);

    old_iter_pressure_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(mesh, MixedPtr<FE_P_disc>(0));
    old_iter_pressure.set(old_iter_pressure_ptr_, 0.0);

    old_div_u_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(eq_data.mechanics_->eq_fields().output_div_ptr->get_dofhandler());
    old_div_u.set(old_div_u_ptr_, 0.0);

    flow_source.set(Model<3, FieldValue<3>::Scalar>::create(
        fn_fluid_source(time_),
        alpha,
        beta,
        eq_data.flow_->eq_fields().field_edge_pressure,
        old_pressure,
        eq_data.mechanics_->eq_fields().output_divergence,
        old_div_u
        ), 0.0);
}

                                    

HM_Iterative::HM_Iterative(Mesh &mesh, Input::Record in_record)
: DarcyFlowInterface(mesh, in_record),
  IterativeCoupling(in_record),
  flow_potential_assembly_(nullptr),
  residual_assembly_(nullptr)
{
	START_TIMER("HM constructor");
    using namespace Input;

    time_ = new TimeGovernor(in_record.val<Record>("time"));
    ASSERT( time_->is_default() == false ).error("Missing key 'time' in Coupling_Iterative.");
    
    // setup flow equation
    Record flow_rec = in_record.val<Record>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    eq_data_.flow_ = std::make_shared<DarcyLMH>(*mesh_, flow_rec, time_);
    
    // setup mechanics
    Record mech_rec = in_record.val<Record>("mechanics_equation");
    eq_data_.mechanics_ = std::make_shared<Elasticity>(*mesh_, mech_rec, this->time_);
    eq_data_.mechanics_->initialize();
    
    // setup coupling fields and finish initialization of flow
    eq_data_.mechanics_->eq_fields()["cross_section"].copy_from(eq_data_.flow_->eq_fields()["cross_section"]);
    eq_data_.flow_->eq_fields() += eq_data_.mechanics_->eq_fields()["cross_section_updated"];
    eq_data_.flow_->eq_fields() += eq_data_.mechanics_->eq_fields()["stress"];
    eq_data_.flow_->eq_fields() += eq_data_.mechanics_->eq_fields()["von_mises_stress"];
    eq_data_.flow_->eq_fields() += eq_data_.mechanics_->eq_fields()["mean_stress"];
    eq_data_.flow_->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "flow") )
        WarningOut() << ss.str();

    eq_fields_ += *eq_data_.flow_->eq_fields().field("pressure_edge");

    this->eq_fieldset_ = &eq_fields_;
    
    // setup input fields
    eq_fields_.set_input_list( in_record.val<Input::Array>("input_fields"), time() );

    eq_fields_.initialize(*mesh_, eq_data_, time_, input_record_.val<double>("iteration_parameter"));
    eq_data_.mechanics_->set_potential_load(eq_fields_.pressure_potential, eq_fields_.ref_pressure_potential);

    eq_fields_.add_coords_field();
}


void HM_Iterative::initialize()
{
    flow_potential_assembly_ = new GenericAssembly<FlowPotentialAssemblyHM>(&eq_fields_, &eq_data_);
    residual_assembly_ = new GenericAssembly<ResidualAssemblyHM>(&eq_fields_, &eq_data_);

    Input::Array user_fields_arr;
    if (input_record_.opt_val("user_fields", user_fields_arr)) {
        FieldSet sham_eq_output; // only for correct call of init_user_fields method
       	this->init_user_fields(user_fields_arr, time().step().end(), sham_eq_output);
    }
}


template<int dim, class Value>
void copy_field(const FieldCommon &from_field_common, FieldFE<dim, Value> &to_field)
{
    auto dh = to_field.get_dofhandler();
    auto vec = to_field.vec();
    Field<dim,Value> from_field;
    from_field.copy_from(from_field_common);
    
    for ( auto cell : dh->own_range() )
        vec.set( cell.local_idx(), from_field.value(cell.elm().centre(), cell.elm()) );
}



void HM_Iterative::zero_time_step()
{
    eq_fields_.set_time(time_->step(), LimitSide::right);
    std::stringstream ss;
    if ( FieldCommon::print_message_table(ss, "coupling_iterative") )
        WarningOut() << ss.str();
    
    eq_data_.mechanics_->update_output_fields(); // init field values for use in flow
    eq_data_.flow_->zero_time_step();
    update_potential();
    eq_data_.mechanics_->zero_time_step();
    
    copy_field(*eq_data_.flow_->eq_fields().field("pressure_p0"), *eq_fields_.old_iter_pressure_ptr_);
    copy_field(eq_data_.mechanics_->eq_fields().output_divergence, *eq_fields_.old_div_u_ptr_);
    eq_fields_.old_iter_pressure.set_time_result_changed();
}


void HM_Iterative::update_solution()
{
    time_->next_time();
    time_->view("HM");
    eq_fields_.set_time(time_->step(), LimitSide::right);

    solve_step();
}

void HM_Iterative::solve_iteration()
{
    // pass displacement (divergence) to flow
    // and solve flow problem
    update_flow_fields();
    eq_data_.flow_->solve_time_step(false);
    
    // pass pressure to mechanics and solve mechanics
    update_potential();
    eq_data_.mechanics_->solve_linear_system();
}


void HM_Iterative::update_after_iteration()
{
    eq_data_.mechanics_->update_output_fields();
    copy_field(*eq_data_.flow_->eq_fields().field("pressure_p0"), *eq_fields_.old_iter_pressure_ptr_);
    eq_fields_.old_iter_pressure.set_time_result_changed();
}


void HM_Iterative::update_after_converged()
{
    eq_data_.flow_->accept_time_step();
    eq_data_.flow_->output_data();
    eq_data_.mechanics_->output_data();
    
    copy_field(eq_data_.mechanics_->eq_fields().output_divergence, *eq_fields_.old_div_u_ptr_);
}


void HM_Iterative::update_potential()
{
    auto ref_potential_vec_ = eq_fields_.ref_potential_ptr_->vec();
    auto dh = eq_fields_.ref_potential_ptr_->get_dofhandler();

    ref_potential_vec_.zero_entries();
    flow_potential_assembly_->assemble(dh);
    
    ref_potential_vec_.local_to_ghost_begin();
    ref_potential_vec_.local_to_ghost_end();
    eq_fields_.pressure_potential.set_time_result_changed();
    eq_fields_.ref_pressure_potential.set_time_result_changed();
    eq_data_.mechanics_->set_potential_load(eq_fields_.pressure_potential, eq_fields_.ref_pressure_potential);
}


void HM_Iterative::update_flow_fields()
{
    eq_fields_.beta.set_time_result_changed();
    eq_fields_.flow_source.set_time_result_changed();
    eq_data_.flow_->set_extra_storativity(eq_fields_.beta);
    eq_data_.flow_->set_extra_source(eq_fields_.flow_source);
}


void HM_Iterative::compute_iteration_error(double& abs_error, double& rel_error)
{
    auto dh = eq_fields_.old_iter_pressure_ptr_->get_dofhandler();
    eq_data_.p_dif2 = 0;
    eq_data_.p_norm2 = 0;

    residual_assembly_->assemble(dh);
    
    double send_data[] = { eq_data_.p_dif2, eq_data_.p_norm2 };
    double recv_data[2];
    MPI_Allreduce(&send_data, &recv_data, 2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    abs_error = sqrt(recv_data[0]);
    rel_error = abs_error / (sqrt(recv_data[1]) + std::numeric_limits<double>::min());
    
    MessageOut().fmt("HM Iteration {} abs. difference: {}  rel. difference: {}\n"
                         "--------------------------------------------------------",
                         iteration(), abs_error, rel_error);

    if(iteration() >= max_it_ && (abs_error > a_tol_ || rel_error > r_tol_))
        THROW(ExcSolverDiverge() << EI_Reason("Reached max_it."));
}



HM_Iterative::~HM_Iterative() {
	eq_data_.flow_.reset();
    eq_data_.mechanics_.reset();
    if (flow_potential_assembly_ != nullptr) delete flow_potential_assembly_;
    if (residual_assembly_ != nullptr) delete residual_assembly_;
}



