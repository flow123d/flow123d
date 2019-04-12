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
#include "fields/field_fe.hh"         // for create_field()


FLOW123D_FORCE_LINK_IN_CHILD(coupling_iterative);


namespace it = Input::Type;



const it::Record & HM_Iterative::get_input_type() {
    return it::Record("Coupling_Iterative",
            "Record with data for iterative coupling of flow and mechanics.\n")
        .derive_from( DarcyFlowInterface::get_input_type() )
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
        .declare_key( "a_tol", it::Double(0), it::Default("1e-7"),
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
    
    *this += pressure_potential.name("pressure_potential")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_result);
}

                                    

HM_Iterative::HM_Iterative(Mesh &mesh, Input::Record in_record)
: DarcyFlowInterface(mesh, in_record)
{
	START_TIMER("HM constructor");
    using namespace Input;
    
    // setup flow equation
    Record flow_rec = in_record.val<Record>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    flow_ = std::make_shared<RichardsLMH>(*mesh_, flow_rec);
    flow_->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "flow") )
        WarningOut() << ss.str();
    
    // setup mechanics
    Record mech_rec = in_record.val<Record>("mechanics_equation");
    mechanics_ = std::make_shared<Elasticity>(*mesh_, mech_rec);
    mechanics_->data()["cross_section"].copy_from(flow_->data()["cross_section"]);
    mechanics_->initialize();
    
    // read parameters controlling the iteration
    beta_ = in_record.val<double>("iteration_parameter");
    min_it_ = in_record.val<unsigned int>("min_it");
    max_it_ = in_record.val<unsigned int>("max_it");
    a_tol_ = in_record.val<double>("a_tol");
    r_tol_ = in_record.val<double>("r_tol");

    this->eq_data_ = &data_;
    
    this->time_ = &flow_->time();
    
    // synchronize time marks of flow and mechanics
    for (auto mark = TimeGovernor::marks().begin(flow_->mark_type()); mark != TimeGovernor::marks().end(flow_->mark_type()); ++mark )
        TimeGovernor::marks().add( TimeMark(mark->time(), mechanics_->time().equation_fixed_mark_type()) );
    for (auto mark = TimeGovernor::marks().begin(mechanics_->mark_type()); mark != TimeGovernor::marks().end(mechanics_->mark_type()); ++mark )
        TimeGovernor::marks().add( TimeMark(mark->time(), flow_->time().equation_fixed_mark_type()) );
    
    data_.set_mesh(*mesh_);
    
    // setup input fields
    data_.set_input_list( in_record.val<Input::Array>("input_fields"), time() );
    
    potential_vec_.resize(mesh_->n_elements());
    potential_ptr_ = create_field<3, FieldValue<3>::Scalar>(potential_vec_, *mesh_, 1);
    data_.pressure_potential.set_field(mesh_->region_db().get_region_set("ALL"), potential_ptr_);
    mechanics_->set_potential_load(data_.pressure_potential);
}


void HM_Iterative::initialize()
{
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
}


void HM_Iterative::update_solution()
{
    unsigned it = 0;
    double difference = 0;
    double init_difference = 1;
    
    data_.set_time(time_->step(), LimitSide::right);
    mechanics_->next_time();
    
    while (it < min_it_ || 
           (it < max_it_ && difference > a_tol_ && difference/init_difference > r_tol_)
          )
    {
        it++;
        
        flow_->update_solution();
        // TODO: pass pressure to mechanics
        update_potential();
        
        mechanics_->solve_linear_system();
        // TODO: pass displacement (divergence) to flow
        // TODO: compute difference of iterates
    }
    mechanics_->output_data();
}


void HM_Iterative::update_potential()
{
    double difference2 = 0, norm2 = 0;
    for (auto ele : mesh_->elements_range())
    {
        double alpha = data_.alpha.value(ele.centre(), ele);
        double pressure = flow_->output_fields().field_ele_pressure.value(ele.centre(), ele);
        double potential = -alpha*pressure;
        difference2 += pow(potential_vec_[ele.idx()] - potential,2);
        norm2 += pow(potential,2);
        potential_vec_[ele.idx()] = potential;
    }
    double dif2norm;
    if (norm2 == 0)
        dif2norm = (difference2 == 0)?0:numeric_limits<double>::max();
    else 
        dif2norm = sqrt(difference2/norm2);
    fill_output_data(potential_vec_, potential_ptr_);
    DebugOut() << "Relative potential difference: " << dif2norm << endl;
    if (dif2norm > numeric_limits<double>::epsilon())
    {
        data_.pressure_potential.set_time_result_changed();
        mechanics_->set_potential_load(data_.pressure_potential);
    }
}


const MH_DofHandler & HM_Iterative::get_mh_dofhandler()
{ 
    return flow_->get_mh_dofhandler(); 
}


HM_Iterative::~HM_Iterative() {
	flow_.reset();
    mechanics_.reset();
}



