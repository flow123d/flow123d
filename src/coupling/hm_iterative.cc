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
#include "flow/richards_lmh.hh"
#include "fields/field_fe.hh"         // for create_field_fe()


FLOW123D_FORCE_LINK_IN_CHILD(coupling_iterative)


namespace it = Input::Type;


const it::Record & HM_Iterative::get_input_type() {
    return it::Record("Coupling_Iterative",
            "Record with data for iterative coupling of flow and mechanics.\n")
        .derive_from( DarcyFlowInterface::get_input_type() )
        .copy_keys(EquationBase::record_template())
        .copy_keys(IterativeCoupling::record_template())
		.declare_key("flow_equation", RichardsLMH::get_input_type(),
		        it::Default::obligatory(),
				"Flow equation, provides the velocity field as a result.")
		.declare_key("mechanics_equation", Elasticity::get_input_type(),
				"Mechanics, provides the displacement field.")
        .declare_key("input_fields", it::Array(
		        HM_Iterative::EqData()
		            .make_field_descriptor_type("Coupling_Iterative")),
		        IT::Default::obligatory(),
		        "Input fields of the HM coupling.")
        .declare_key( "iteration_parameter", it::Double(), it::Default("1"),
                "Tuning parameter for iterative splitting. Its default value"
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


HM_Iterative::EqData::EqData()
{
    *this += alpha.name("biot_alpha")
                     .units(UnitSI().dimensionless())
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += density.name("fluid_density")
                     .units(UnitSI().kg().m(-3))
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += gravity.name("gravity")
                     .units(UnitSI().m().s(-2))
                     .input_default("9.81")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += beta.name("relaxation_beta")
                     .units(UnitSI().dimensionless())
                     .flags(FieldFlag::equation_external_output);
    
    *this += pressure_potential.name("pressure_potential")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_result);
    
    *this += flow_source.name("extra_flow_source")
                     .units(UnitSI().s(-1))
                     .flags(FieldFlag::equation_result);
}


void HM_Iterative::EqData::initialize(Mesh &mesh)
{
    // initialize coupling fields with FieldFE
    set_mesh(mesh);
    
    potential_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(mesh, MixedPtr<FE_CR>());
    pressure_potential.set_field(mesh.region_db().get_region_set("ALL"), potential_ptr_);
    
    beta_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(mesh, MixedPtr<FE_P_disc>(0));
    beta.set_field(mesh.region_db().get_region_set("ALL"), beta_ptr_);
    
    flow_source_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(beta_ptr_->get_dofhandler());
    flow_source.set_field(mesh.region_db().get_region_set("ALL"), flow_source_ptr_);
    
    old_pressure_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(beta_ptr_->get_dofhandler());
    old_iter_pressure_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(beta_ptr_->get_dofhandler());
    div_u_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(beta_ptr_->get_dofhandler());
    old_div_u_ptr_ = create_field_fe<3, FieldValue<3>::Scalar>(beta_ptr_->get_dofhandler());
}

                                    

HM_Iterative::HM_Iterative(Mesh &mesh, Input::Record in_record)
: DarcyFlowInterface(mesh, in_record),
  IterativeCoupling(in_record)
{
	START_TIMER("HM constructor");
    using namespace Input;

    time_ = new TimeGovernor(in_record.val<Record>("time"));
    ASSERT( time_->is_default() == false ).error("Missing key 'time' in Coupling_Iterative.");
    
    // setup flow equation
    Record flow_rec = in_record.val<Record>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    flow_ = std::make_shared<RichardsLMH>(*mesh_, flow_rec, time_);
    flow_->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "flow") )
        WarningOut() << ss.str();
    
    // setup mechanics
    Record mech_rec = in_record.val<Record>("mechanics_equation");
    mechanics_ = std::make_shared<Elasticity>(*mesh_, mech_rec, this->time_);
    mechanics_->data()["cross_section"].copy_from(flow_->data()["cross_section"]);
    mechanics_->initialize();
    
    // read parameters controlling the iteration
    beta_ = in_record.val<double>("iteration_parameter");

    this->eq_data_ = &data_;
    
    // setup input fields
    data_.set_input_list( in_record.val<Input::Array>("input_fields"), time() );

    data_.initialize(*mesh_);
    mechanics_->set_potential_load(data_.pressure_potential);
}


void HM_Iterative::initialize()
{
}


template<int dim, class Value>
void copy_field(const FieldCommon &from_field_common, FieldFE<dim, Value> &to_field)
{
    auto dh = to_field.get_dofhandler();
    auto vec = to_field.get_data_vec();
    Field<dim,Value> from_field;
    from_field.copy_from(from_field_common);
    
    for ( auto cell : dh->own_range() )
        vec[cell.local_idx()] = from_field.value(cell.elm().centre(), cell.elm());
}



void HM_Iterative::zero_time_step()
{
    data_.set_time(time_->step(), LimitSide::right);
    std::stringstream ss;
    if ( FieldCommon::print_message_table(ss, "coupling_iterative") )
        WarningOut() << ss.str();
    
    flow_->zero_time_step();
    update_potential();
    mechanics_->zero_time_step();
    
    copy_field(*flow_->data().field("pressure_p0"), *data_.old_pressure_ptr_);
    copy_field(*flow_->data().field("pressure_p0"), *data_.old_iter_pressure_ptr_);
    copy_field(mechanics_->data().output_divergence, *data_.div_u_ptr_);
}


void HM_Iterative::update_solution()
{
    time_->next_time();
    time_->view("HM");
    data_.set_time(time_->step(), LimitSide::right);

    solve_step();
}

void HM_Iterative::solve_iteration()
{
    // pass displacement (divergence) to flow
    // and solve flow problem
    update_flow_fields();
    flow_->solve_time_step(false);
    
    // pass pressure to mechanics and solve mechanics
    update_potential();
    mechanics_->solve_linear_system();
}


void HM_Iterative::update_after_iteration()
{
    mechanics_->update_output_fields();
    copy_field(mechanics_->data().output_divergence, *data_.div_u_ptr_);
    copy_field(*flow_->data().field("pressure_p0"), *data_.old_iter_pressure_ptr_);
}


void HM_Iterative::update_after_converged()
{
    flow_->accept_time_step();
    flow_->output_data();
    mechanics_->output_data();
    
    copy_field(*flow_->data().field("pressure_p0"), *data_.old_pressure_ptr_);
    copy_field(mechanics_->data().output_divergence, *data_.old_div_u_ptr_);
}


void HM_Iterative::update_potential()
{
    auto potential_vec_ = data_.potential_ptr_->get_data_vec();
    auto dh = data_.potential_ptr_->get_dofhandler();
    Field<3, FieldValue<3>::Scalar> field_edge_pressure;
    field_edge_pressure.copy_from(*flow_->data().field("pressure_edge"));
    for ( auto ele : dh->local_range() )
    {
        auto elm = ele.elm();
        LocDofVec dof_indices = ele.get_loc_dof_indices();
        for ( auto side : ele.side_range() )
        {
            double alpha = data_.alpha.value(side.centre(), elm);
            double density = data_.density.value(side.centre(), elm);
            double gravity = data_.gravity.value(side.centre(), elm);
            double pressure = field_edge_pressure.value(side.centre(), elm);
            double potential = -alpha*density*gravity*pressure;
        
            potential_vec_[dof_indices[side.side_idx()]] = potential;
        }
    }
    
    data_.pressure_potential.set_time_result_changed();
    mechanics_->set_potential_load(data_.pressure_potential);
}


void HM_Iterative::update_flow_fields()
{
    auto beta_vec = data_.beta_ptr_->get_data_vec();
    auto src_vec = data_.flow_source_ptr_->get_data_vec();
    auto dh = data_.beta_ptr_->get_dofhandler();
    Field<3,FieldValue<3>::Scalar> field_ele_pressure;
    field_ele_pressure.copy_from(*flow_->data().field("pressure_p0"));
    for ( auto ele : dh->local_range() )
    {
        auto elm = ele.elm();
        
        double alpha = data_.alpha.value(elm.centre(), elm);
        double young = mechanics_->data().young_modulus.value(elm.centre(), elm);
        double poisson = mechanics_->data().poisson_ratio.value(elm.centre(), elm);
        double beta = beta_ * 0.5*alpha*alpha/(2*lame_mu(young, poisson)/elm.dim() + lame_lambda(young, poisson));
        
        double old_p = data_.old_pressure_ptr_->value(elm.centre(), elm);
        double p = field_ele_pressure.value(elm.centre(), elm);
        double div_u = data_.div_u_ptr_->value(elm.centre(), elm);
        double old_div_u = data_.old_div_u_ptr_->value(elm.centre(), elm);
        double src = (beta*(p-old_p) + alpha*(old_div_u - div_u)) / time_->dt();
        
        beta_vec[ele.local_idx()] = beta;
        src_vec[ele.local_idx()] = src;
    }
    
    data_.beta.set_time_result_changed();
    data_.flow_source.set_time_result_changed();
    flow_->set_extra_storativity(data_.beta);
    flow_->set_extra_source(data_.flow_source);
}


void HM_Iterative::compute_iteration_error(double& abs_error, double& rel_error)
{
    auto dh = data_.beta_ptr_->get_dofhandler();
    double p_dif2 = 0, p_norm2 = 0;
    Field<3,FieldValue<3>::Scalar> field_ele_pressure;
    field_ele_pressure.copy_from(*flow_->data().field("pressure_p0"));
    for (auto cell : dh->own_range())
    {
        auto elm = cell.elm();
        double new_p = field_ele_pressure.value(elm.centre(), elm);
        double old_p = data_.old_iter_pressure_ptr_->value(elm.centre(), elm);
        p_dif2 += pow(new_p - old_p, 2)*elm.measure();
        p_norm2 += pow(old_p, 2)*elm.measure();
    }
    
    double send_data[] = { p_dif2, p_norm2 };
    double recv_data[2];
    MPI_Allreduce(&send_data, &recv_data, 2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    abs_error = sqrt(recv_data[0]);
    rel_error = abs_error / sqrt(recv_data[1]);
    
    MessageOut().fmt("HM Iteration {} abs. difference: {}  rel. difference: {}\n"
                         "--------------------------------------------------------",
                         iteration(), abs_error, rel_error);
}



double HM_Iterative::last_t()
{
    return flow_->last_t();
}


HM_Iterative::~HM_Iterative() {
	flow_.reset();
    mechanics_.reset();
}



