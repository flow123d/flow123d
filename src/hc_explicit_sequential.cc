/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief
 *
 *  @author Jan Brezina
 */

#include "hc_explicit_sequential.hh"
#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "transport/transport_operator_splitting.hh"
#include "transport/transport.h"
#include "transport/transport_dg.hh"
#include "equation.hh"
#include "time_marks.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "io/output.h"
#include "main.h"

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR UNSTEADY SATURATED FLOW
 */
HC_ExplicitSequential::HC_ExplicitSequential(ProblemType problem_type)
{
    type_=problem_type;

    // Initialize Time Marks
    main_time_marks = new TimeMarks();

    // Material Database
    const string& material_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Material", NULL ));
    material_database = new MaterialDatabase(material_file_name);

    // Read mesh
    mesh = new Mesh();
    const string& mesh_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr("Input", "Mesh", NULL));
    MeshReader* meshReader = new GmshMeshReader();
    meshReader->read(mesh_file_name, mesh);
    mesh->setup_topology();
    mesh->setup_materials(*material_database);
    Profiler::instance()->set_task_size(mesh->n_elements());

    // setup water flow object
    switch (problem_type) {
        case STEADY_SATURATED:
            water=new DarcyFlowMH_Steady(*main_time_marks, *mesh, *material_database);
            break;
        case UNSTEADY_SATURATED:
            DBGMSG("Unsteady\n");
            water=new DarcyFlowMH_Unsteady(*main_time_marks, *mesh, *material_database);
            break;
        case UNSTEADY_SATURATED_LMH:
            water=new DarcyFlowLMH_Unsteady(*main_time_marks, *mesh, *material_database);
            break;
        default:
            xprintf(Err,"Wrong problem type: %d",problem_type);
            break;
    }
    // object for water postprocessing and output

    water_output = new DarcyFlowMHOutput(water);

    // optionally setup transport objects
    if ( OptGetBool("Transport", "Transport_on", "no") ) {
        char *transport_type = OptGetStr("Transport", "Transport_type", "explicit");
        if (strcmp(transport_type, "explicit") == 0)
        {
            transport_reaction = new TransportOperatorSplitting(*main_time_marks, *mesh, *material_database);
        }
        else if (strcmp(transport_type, "implicit") == 0)
        {
            transport_reaction = new TransportDG(*main_time_marks, *mesh, *material_database);
        }
        else
        {
            xprintf(PrgErr,"Value of parameter: [Transport] Transport_type is neither \"explicit\" nor \"implicit\".\n");
        }
    } else {
        transport_reaction = new TransportNothing(*main_time_marks, *mesh, *material_database);
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

    // following should be specified in constructor:
    // value for velocity interpolation :
    // theta = 0     velocity from beginning of transport interval (fully explicit method)
    // theta = 0.5   velocity from center of transport interval ( mimic Crank-Nicholson)
    // theta = 1.0   velocity from end of transport interval (partialy explicit scheme)
    const double theta=0.5;

    double velocity_interpolation_time;
    bool velocity_changed;
    Vec velocity_field;


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

    // output initial condition
    water_output->postprocess();
    water_output->output();


    while (! (water->time().is_end() && transport_reaction->time().is_end() ) ) {

        transport_reaction->set_time_step_constrain(water->time().dt());
        // in future here could be re-estimation of transport planed time according to
        // evolution of the velocity field. Consider the case w_dt << t_dt and velocity almost constant in time
        // which suddenly rise in time 3*w_dt. First we the planed transport time step t_dt could be quite big, but
        // in time 3*w_dt we can reconsider value of t_dt to better capture changing velocity.
        velocity_interpolation_time= theta * transport_reaction->planned_time() + (1-theta) * transport_reaction->solved_time();

        // if transport is off, transport should return infinity solved and planned times so that
        // only water branch takes the place
        if (water->solved_time() < velocity_interpolation_time) {
            // solve water over the nearest transport interval
            water->update_solution();
            water_output->postprocess();
            // here possibly save solution from water for interpolation in time

            water_output->output();

            water->choose_next_time();

            velocity_changed = true;
        } else {
            // having information about velocity field we can perform transport step

            // here should be interpolation of the velocity at least if the interpolation time
            // is not close to the solved_time of the water module
            // for simplicity we use only last velocity field
            if (velocity_changed) {
                //DBGMSG("velocity update\n");
                water->get_velocity_seq_vector(velocity_field);
                transport_reaction->set_velocity_field(velocity_field);
                velocity_changed = false;
            }
            transport_reaction->update_solution();
            transport_reaction->output_data();
        }

    }
    xprintf(Msg, "End of simulation at time: %f\n", transport_reaction->solved_time());
}


HC_ExplicitSequential::~HC_ExplicitSequential() {

    delete mesh;
    delete material_database;
    delete main_time_marks;
    delete water;
    delete water_output;
    delete transport_reaction;
}



//-----------------------------------------------------------------------------
// vim: set cindent:
//-----------------------------------------------------------------------------
