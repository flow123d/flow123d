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
 * @file    hc_explicit_sequential.cc
 * @brief   
 * @author  Jan Brezina
 */

#include "hc_explicit_sequential.hh"
#include "flow/darcy_flow_interface.hh"
#include "flow/darcy_flow_mh.hh"
// TODO:
// After having general default values:
// make TransportNoting default for AdvectionProcessBase abstract
// use default "{}" for secondary equation.
// Then we can remove following include.

#include "transport/transport_operator_splitting.hh"
#include "fields/field_common.hh"
#include "transport/heat_model.hh"

#include "fields/field_set.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "system/sys_profiler.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"


FLOW123D_FORCE_LINK_IN_PARENT(transportOperatorSplitting)
FLOW123D_FORCE_LINK_IN_PARENT(concentrationTransportModel)
FLOW123D_FORCE_LINK_IN_PARENT(convectionTransport)
FLOW123D_FORCE_LINK_IN_PARENT(heatModel)

FLOW123D_FORCE_LINK_IN_PARENT(darcy_flow_mh)
FLOW123D_FORCE_LINK_IN_PARENT(darcy_flow_lmh)
FLOW123D_FORCE_LINK_IN_PARENT(richards_lmh)
FLOW123D_FORCE_LINK_IN_PARENT(coupling_iterative)


namespace it = Input::Type;

it::Abstract & CouplingBase::get_input_type() {
	return it::Abstract("Coupling_Base", "The root record of description of particular the problem to solve.")
		.close();
}


const it::Record & HC_ExplicitSequential::get_input_type() {
    return it::Record("Coupling_Sequential",
            "Record with data for a general sequential coupling.\n")
		.derive_from( CouplingBase::get_input_type() )
		.declare_key("description",it::String(),
				"Short description of the solved problem.\n"
				"Is displayed in the main log, and possibly in other text output files.")
		.declare_key("mesh", Mesh::get_input_type(), it::Default::obligatory(),
				"Computational mesh common to all equations.")
		.declare_key("time", TimeGovernor::get_input_type(), it::Default::optional(),
				"Simulation time frame and time step.")
		.declare_key("flow_equation", DarcyFlowInterface::get_input_type(),
		        it::Default::obligatory(),
				"Flow equation, provides the velocity field as a result.")
		.declare_key("solute_equation", AdvectionProcessBase::get_input_type(),
				"Transport of soluted substances, depends on the velocity field from a Flow equation.")
		.declare_key("heat_equation", AdvectionProcessBase::get_input_type(),
		        "Heat transfer, depends on the velocity field from a Flow equation.")
		.close();
}


const int HC_ExplicitSequential::registrar = HC_ExplicitSequential::get_input_type().size();



std::shared_ptr<AdvectionProcessBase> HC_ExplicitSequential::make_advection_process(string process_key)
{
    using namespace Input;
    // setup heat object
    Iterator<AbstractRecord> it = in_record_.find<AbstractRecord>(process_key);

    if (it) {
        auto process = (*it).factory< AdvectionProcessBase, Mesh &, const Input::Record >(*mesh, *it);

        // setup fields
        process->data()["cross_section"]
                .copy_from(water->data()["cross_section"]);
        /*
        if (water_content_saturated_) // only for unsteady Richards water model
            process->data()["porosity"].copy_from(*water_content_saturated_);

        if (water_content_p0_)
            process->data()["water_content"].copy_from(*water_content_p0_);
        else {

        }*/

        FieldCommon *porosity = process->data().field("porosity");
        process->data()["water_content"].copy_from( *porosity );


        process->initialize();
        return process;
    } else {
        return std::make_shared<TransportNothing>(*mesh);
    }
}

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR UNSTEADY SATURATED FLOW
 */
HC_ExplicitSequential::HC_ExplicitSequential(Input::Record in_record)
: in_record_(in_record)
{
	START_TIMER("HC constructor");
    using namespace Input;

    // Read mesh
    {
        START_TIMER("HC read mesh");

   		mesh = BaseMeshReader::mesh_factory( in_record.val<Record>("mesh") );
        
        //getting description for the Profiler
        string description;
        in_record.opt_val<string>("description", description);
         
        Profiler::instance()->set_task_info(
            description,
            //"Description has to be set in main. by different method.",
            mesh->n_elements());
    }

    // setup primary equation - water flow object
    AbstractRecord prim_eq = in_record.val<AbstractRecord>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    water = prim_eq.factory< DarcyFlowInterface, Mesh &, const Input::Record>(*mesh, prim_eq);
    water->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "HC explicit sequential") ) {
        WarningOut() << ss.str();
    }

    RegionSet bulk_set = mesh->region_db().get_region_set("BULK");
    water_content_saturated_ = water->data().field("water_content_saturated");
    if (water_content_saturated_ && water_content_saturated_->field_result( bulk_set ) == result_zeros )
        water_content_saturated_ = nullptr;

    water_content_p0_ = water->data().field("water_content_p0");
    if (water_content_p0_ && water_content_p0_->field_result( bulk_set ) == result_zeros )
        water_content_p0_ = nullptr;

    processes_.push_back(AdvectionData(make_advection_process("solute_equation")));
    processes_.push_back(AdvectionData(make_advection_process("heat_equation")));
}

void HC_ExplicitSequential::advection_process_step(AdvectionData &pdata)
{
    if (pdata.process->time().is_end()) return;

    is_end_all_=false;
    if ( pdata.process->time().step().le(pdata.velocity_time) ) {
        // having information about velocity field we can perform transport step

        // here should be interpolation of the velocity at least if the interpolation time
        // is not close to the solved_time of the water module
        // for simplicity we use only last velocity field
        if (pdata.velocity_changed) {
            //DBGMSG("velocity update\n");
//             std::dynamic_pointer_cast<DarcyMH>(water)->get_velocity_field()->local_to_ghost_data_scatter_begin();
//             std::dynamic_pointer_cast<DarcyMH>(water)->get_velocity_field()->local_to_ghost_data_scatter_end();
//             pdata.process->set_velocity_field( std::dynamic_pointer_cast<DarcyMH>(water)->get_velocity_field() );
            water->get_velocity_field()->local_to_ghost_data_scatter_begin();
            water->get_velocity_field()->local_to_ghost_data_scatter_end();
            pdata.process->set_velocity_field( water->get_velocity_field() );
            pdata.velocity_changed = false;
        }
        if (pdata.process->time().tlevel() == 0) pdata.process->zero_time_step();

        pdata.process->update_solution();
        pdata.process->output_data();
   }
}


/**
 * TODO:
 * - have support for steady problems in TimeGovernor, make Noting problems steady
 * - apply splitting of compute_one_step to particular models
 * - how to set output time marks for steady problems (we need solved time == infinity) but
 *   add no time marks
 * - allow create steady time governor without time marks (at least in nothing models)
 * - pass refference to time marks in time governor constructor?
 */

void HC_ExplicitSequential::run_simulation()
{
    START_TIMER("HC run simulation");
    // following should be specified in constructor:
    // value for velocity interpolation :
    // theta = 0     velocity from beginning of transport interval (fully explicit method)
    // theta = 0.5   velocity from center of transport interval ( mimic Crank-Nicholson)
    // theta = 1.0   velocity from end of transport interval (partialy explicit scheme)
    const double theta=0.5;

    {
        START_TIMER("HC water zero time step");
        water->zero_time_step();
    }


    // following cycle is designed to support independent time stepping of
    // both processes. The question is which value of the water field use to compute a transport step.
    // Meaningful cases are
    //      1) beginning (fully explicit method)
    //      2) center ( mimic Crank-Nicholson)
    //      3) end of the interval (partialy explicit scheme)
    // However with current implementation of the explicit transport on have to assembly transport matrix for
    // every new value of the velocity field. So we have to keep same velocity field over some time interval t_dt
    // which is further split into shorter time intervals ts_dt dictated by the CFL condition.
    // One can consider t_dt as the transport time step and apply one of the previous three cases.
    //
    // The question is how to choose intervals t_dt. That should depend on variability of the velocity field in time.
    // Currently we simply use t_dt == w_dt.

    is_end_all_=false;
    while (! is_end_all_) {
        is_end_all_ = true;

        double water_dt=water->time().estimate_dt();
        if (water->time().is_end()) water_dt = TimeGovernor::inf_time;

        // in future here could be re-estimation of transport planed time according to
        // evolution of the velocity field. Consider the case w_dt << t_dt and velocity almost constant in time
        // which suddenly rise in time 3*w_dt. First we the planed transport time step t_dt could be quite big, but
        // in time 3*w_dt we can reconsider value of t_dt to better capture changing velocity.
        min_velocity_time = TimeGovernor::inf_time;
        for(auto &pdata : processes_) {
            pdata.process->set_time_upper_constraint(water_dt, "Flow time step");
            pdata.velocity_time = theta * pdata.process->planned_time() + (1-theta) * pdata.process->solved_time();
            min_velocity_time = min(min_velocity_time, pdata.velocity_time);
        }

        // printing water and transport times every step
        //MessageOut().fmt("HC_EXPL_SEQ: velocity_interpolation_time: {}, water_time: {} transport time: {}\n",
        //        velocity_interpolation_time, water->time().t(), transport_reaction->time().t());
         
        // if transport is off, transport should return infinity solved and planned times so that
        // only water branch takes the place

        if (! water->time().is_end() ) {
            is_end_all_=false;
            if (water->solved_time() < min_velocity_time) {

                // solve water over the nearest transport interval
                water->update_solution();

                // here possibly save solution from water for interpolation in time

                //water->time().view("WATER");     //show water time governor

                //water->output_data();
                water->choose_next_time();

                for(auto &process : processes_) process.velocity_changed = true;

            }
        }
        advection_process_step(processes_[0]); // solute
        advection_process_step(processes_[1]); // heat
    }
    //MessageOut().fmt("End of simulation at time: {}\n", max(solute->solved_time(), heat->solved_time()));
}


HC_ExplicitSequential::~HC_ExplicitSequential() {
	water.reset();
	for(auto &pdata : processes_) pdata.process.reset();
    delete mesh;
}




//-----------------------------------------------------------------------------
// vim: set cindent:
//-----------------------------------------------------------------------------
