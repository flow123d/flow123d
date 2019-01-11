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


FLOW123D_FORCE_LINK_IN_CHILD(coupling_iterative);


namespace it = Input::Type;



const it::Record & HM_Iterative::get_input_type() {
    return it::Record("Coupling_Iterative",
            "Record with data for iterative coupling of flow and mechanics.\n")
        .derive_from( DarcyFlowInterface::get_input_type() )
		.declare_key("flow_equation", RichardsLMH::get_input_type(),
		        it::Default::obligatory(),
				"Flow equation, provides the velocity field as a result.")
// 		.declare_key("mechanics_equation", Mechanics::get_input_type(),
// 				"Mechanics, provides the displacement field.")
        .declare_key( "iteration_parameter", it::Double(), it::Default("1"),
                "Tuning parameter for iterative splitting." )
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



HM_Iterative::HM_Iterative(Mesh &mesh, Input::Record in_record)
: DarcyFlowInterface(mesh, in_record)
{
	START_TIMER("HM constructor");
    using namespace Input;

    // setup flow equation
    Record prim_eq = in_record.val<Record>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    flow_ = std::make_shared<RichardsLMH>(*mesh_, prim_eq);
    flow_->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "HM iterative") ) {
        WarningOut() << ss.str();
    }
    
    // setup mechanics
    // TODO: supply real equation
    mechanics_ = std::make_shared<EquationBase>();
    
    // read parameters controlling the iteration
    beta_ = in_record.val<double>("iteration_parameter");
    min_it_ = in_record.val<unsigned int>("min_it");
    max_it_ = in_record.val<unsigned int>("max_it");
    a_tol_ = in_record.val<double>("a_tol");
    r_tol_ = in_record.val<double>("r_tol");
    
    this->eq_data_ = std::make_shared<FieldSet>().get();
    this->time_ = &flow_->time();
}


void HM_Iterative::zero_time_step()
{
    flow_->zero_time_step();
}


void HM_Iterative::update_solution()
{
    unsigned it = 0;
    double difference = 0;
    double init_difference = 1;
    
    while (it < min_it_ || 
           (it < max_it_ && difference > a_tol_ && difference/init_difference > r_tol_)
          )
    {
        it++;
        
        mechanics_->update_solution();
        // TODO: pass displacement (divergence) to flow
        flow_->update_solution();
        // TODO: pass pressure to mechanics
        // TODO: compute difference of iterates
    }
}


HM_Iterative::~HM_Iterative() {
	flow_.reset();
    mechanics_.reset();
}



