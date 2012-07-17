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
 * @file main.cc
 * @brief This file should contain only creation of Application object.
 *
 */


#include <petsc.h>
#include <sstream>

#include "system/system.hh"
#include "io_namehandler.hh"
#include "hc_explicit_sequential.hh"

#include "main.h"
#include "io/read_ini.h"

#include "rev_num.h"

/// named version of the program
#define _VERSION_   "1.7.0_dev"

static void main_convert_to_output();

/**
 * @brief Main flow initialization
 * @param[in] argc       command line argument count
 * @param[in] argv       command line arguments
 * @param[out] goal      Flow computation goal
 * @param[out] ini_fname Init file name
 *
 * TODO: this parsing function should be in main.cc
 *
 */
void parse_cmd_line(const int argc, char * argv[],  string &ini_fname) {
    const char USAGE_MSG[] = "\
    Wrong program parameters.\n\
    Usage: flow123d [options] ini_file\n\
    Options:\n\
    -s       Compute MH problem (Obsolete)\n\
             Source files have to be in the current directory.\n\
    -S       Compute MH problem\n\
             Source files have to be in the same directory as ini file.\n\
    -i       String used to change the 'variable' ${INPUT} in the file path.\n\
    -o       Absolute path to output directory.\n\
    -l file  Set base name of log files or turn logging off if no name is given.\n";

    xprintf(MsgLog, "Parsing program parameters ...\n");

    // Check command line arguments
    if ((argc >= 3) && (strlen(argv[1]) == 2) && (argv[1][0] == '-')) {
        std::string ini_argument ( argv[2] );
        std::string ini_dir;

        // Try to find absolute or relative path in fname
        int delim_pos=ini_argument.find_last_of(DIR_DELIMITER);
        if (delim_pos < ini_argument.npos) {
            // It seems, that there is some path in fname ... separate it
            ini_dir=ini_argument.substr(0,delim_pos);
            ini_fname=ini_argument.substr(delim_pos+1); // till the end
        } else {
            ini_dir=".";
            ini_fname=ini_argument;
        }

        switch (argv[ 1 ][ 1 ]) {
            case 's':
                ini_fname=ini_argument;
                break;
            case 'S':
                xchdir(ini_dir.c_str());
                break;
            default:
                //xprintf(UsrErr, USAGE_MSG);   // Caused crash of flow123d
                xprintf(UsrErr,"%s", USAGE_MSG);
        }

    }
}

//=============================================================================

/**
 *  FUNCTION "MAIN"
 */
int main(int argc, char **argv) {
    std::string ini_fname;

    F_ENTRY;

    parse_cmd_line(argc, argv,  ini_fname); // command-line parsing

    // temporary moving PETSC stuff here from system.hh
    // should be made better in JB_1.7.input
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

    // determine logfile name or switch it off
    PetscBool flg;
    char file[PETSC_MAX_PATH_LEN];     /* log file name */
    stringstream log_name;

    PetscOptionsGetString(PETSC_NULL,"-l",file,PETSC_MAX_PATH_LEN,&flg);
    if (flg == PETSC_TRUE) {
        if (file[0] == '\n') {
           // -l option without given name -> turn logging off
        } else {
           // given log name
           log_name << string(file) <<  "." << sys_info.my_proc << ".log";

        }
    } else {
        // use default name
        log_name << "flow123."<< sys_info.my_proc << ".log";
    }

    
    system_init(PETSC_COMM_WORLD, IONameHandler::get_instance()->get_output_file_name(log_name.str()) ); 
    OptionsInit(ini_fname.c_str()); // Read options/ini file into database
    system_set_from_options();

    
    Profiler::initialize(MPI_COMM_WORLD, IONameHandler::get_instance()->get_output_file_name(""));

    START_TIMER("WHOLE PROGRAM");

    // Say Hello
    
    // make strings from macros in order to check type
    string version(_VERSION_);
    string revision(REVISION);
    
    xprintf(Msg, "This is FLOW-1-2-3, version %s rev: %s\n", version.c_str(),revision.c_str());
    xprintf(Msg, "Built on %s at %s.\n", __DATE__, __TIME__);

    ProblemType type = (ProblemType) OptGetInt("Global", "Problem_type", NULL);
    switch (type) {
    case CONVERT_TO_OUTPUT:
        main_convert_to_output();
        break;
    case STEADY_SATURATED:
    case UNSTEADY_SATURATED:
    case UNSTEADY_SATURATED_LMH: {
        HC_ExplicitSequential *problem = new HC_ExplicitSequential(type);
        problem->run_simulation();
        delete problem;
        break;
    }
    case PROBLEM_DENSITY:
        // main_compute_mh_density(problem);
        xprintf(UsrErr,"Density driven model not yet reimplemented.");
        break;
    }

    // Say Goodbye
    return xterminate(false);
}

/**
 * FUNCTION "MAIN" FOR CONVERTING FILES TO POS
 */
void main_convert_to_output() {
    // TODO: implement output of input data fields
    // Fields to output:
    // 1) volume data (simple)
    //    sources (Darcy flow and transport), initial condition, material id, partition id
    // 2) boundary data (needs "virtual fractures" in output mesh)
    //    flow and transport bcd

    xprintf(Err, "Not implemented yet in this version\n");
}
#if 0
/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM
 */
void main_compute_mh(struct Problem *problem) {
    int type=OptGetInt("Global", "Problem_type", NULL);
    switch (type) {
        case STEADY_SATURATED:
            main_compute_mh_steady_saturated(problem);
            break;
        case UNSTEADY_SATURATED:
            main_compute_mh_unsteady_saturated(problem);
            break;
        case PROBLEM_DENSITY:
           // main_compute_mh_density(problem);
            break;
        default:
            xprintf(UsrErr,"Unsupported problem type: %d.",type);
    }
}

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR UNSTEADY SATURATED FLOW
 */
void main_compute_mh_unsteady_saturated(struct Problem *problem)
{

    const string& mesh_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr("Input", "Mesh", NULL));
    MeshReader* meshReader = new GmshMeshReader();

    Mesh* mesh = new Mesh();
    meshReader->read(mesh_file_name, mesh);
    mesh->setup_topology();
    mesh->setup_materials(* problem->material_database);
    Profiler::instance()->set_task_size(mesh->n_elements());
    OutputTime *output_time;
    TimeMarks * main_time_marks = new TimeMarks();
    int i;

    // setup output
    string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", "\\"));
    output_time = new OutputTime(mesh, output_file);

    DarcyFlowMH *water = new DarcyFlowLMH_Unsteady(main_time_marks,mesh, problem->material_database);
    DarcyFlowMHOutput *water_output = new DarcyFlowMHOutput(water);
    const TimeGovernor &water_time=water->get_time();

    // set output time marks
    TimeMark::Type output_mark_type = main_time_marks->new_strict_mark_type();
    main_time_marks->add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), water_time.end_time(), output_mark_type );

    while (! water_time.is_end()) {
        water->compute_one_step();
        water_output->postprocess();

        if ( main_time_marks->is_current(water_time, output_mark_type) )  {
            output_time->get_data_from_mesh();
            // call output_time->register_node_data(name, unit, 0, data) to register other data on nodes
            // call output_time->register_elem_data(name, unit, 0, data) to register other data on elements
            output_time->write_data(water_time.t());
            output_time->free_data_from_mesh();
        }
    }
}

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR STEADY SATURATED FLOW
 */
void main_compute_mh_steady_saturated(struct Problem *problem)
{
    const string& mesh_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr("Input", "Mesh", NULL));
    MeshReader* meshReader = new GmshMeshReader();

    Mesh* mesh = new Mesh();
    meshReader->read(mesh_file_name, mesh);
    mesh->setup_topology();
    mesh->setup_materials(* problem->material_database);
    Profiler::instance()->set_task_size(mesh->n_elements());

    TimeMarks * main_time_marks = new TimeMarks();

    /*
       Mesh* mesh;
       ElementIter elm;
       struct Side *sde;
       FILE *out;
       int i;
       mesh=problem->mesh;
     */

    problem->water=new DarcyFlowMH_Steady(main_time_marks, mesh, problem->material_database);
    // Pointer at Output should be in this object
    DarcyFlowMHOutput *water_output = new DarcyFlowMHOutput(problem->water);

    problem->water->compute_one_step();

    if (OptGetBool("Transport", "Transport_on", "no") == true) {
        problem->otransport = new ConvectionTransport(problem->material_database, mesh);
        problem->transport_os = new TransportOperatorSplitting(problem->material_database, mesh);
    }


    water_output->postprocess();

    /* Write static data to output file */
	string out_fname =  IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", NULL));
	Output *output = new Output(mesh, out_fname);
	output->get_data_from_mesh();
	// call output->register_node_data(name, unit, data) here to register other data on nodes
	// call output->register_elem_data(name, unit, data) here to register other data on elements
	output->write_data();
	output->free_data_from_mesh();
	delete output;

    // pracovni vystup nekompatibilniho propojeni
    // melo by to byt ve water*
    /*
       {
            ElementIter ele;
            Element *ele2;
            int ngi;
            Neighbour *ngh;
            Mesh* mesh = problem->mesh;

            double sum1,sum2;
            DarcyFlowMH *w=problem->water;
            double *x = w->schur0->vx;


            FOR_ELEMENTS_IT( ele ) {
                FOR_ELM_NEIGHS_VV( ele, ngi ) {
                    ngh = ele->neigh_vv[ ngi ];
                    // get neigbour element, and set appropriate column
               ele2 = ( ngh->element[ 0 ] == &(*ele) ) ? ngh->element[ 1 ] : ngh->element[ 0 ];

               double out_flux=0.0;
               for(i=0;i<=ele->dim;i++) out_flux+=x[w->side_row_4_id[ele->side[i]->id]];
               xprintf(Msg,"El 1 (%f,%f) %d %f %g\n",
                       ele->centre[0],
                       ele->centre[1],
                       ele->dim,
                       x[w->el_row_4_id[ele->id]],out_flux);

               out_flux=0.0;
               for(i=0;i<=ele2->dim;i++) out_flux+=x[w->side_row_4_id[ele2->side[i]->id]];

               xprintf(Msg,"El 2 (%f,%f) %d %f %g\n",
                                  ele2->centre[0],
                                  ele2->centre[1],
                                  ele2->dim,
                                  x[w->el_row_4_id[ele2->id]],out_flux);
                }
            }
       }

        out = xfopen("pepa.txt","wt");

        FOR_ELEMENTS(elm)
                    elm->aux = 0;

        FOR_ELEMENTS(elm)
            FOR_ELEMENT_SIDES(elm,i)
                            if(elm->side[i]->type == EXTERNAL)
                                    if((elm->side[i]->centre[2] - elm->centre[2]) > 0.0)
                                            if(elm->side[i]->normal[2] != 0.0)
                                                    if((elm->side[i]->centre[2]) > 300.0)
                                                            if((elm->material->id == 2200) || (elm->material->id == 2207) || (elm->material->id == 2212) || (elm->material->id == 2217) || (elm->material->id == 9100) || (elm->material->id == 9107) || (elm->material->id == 9112) || (elm->material->id == 9117))
                                                                    elm->aux = 1;

        FOR_ELEMENTS(elm)
                    if(elm->aux == 1)
                            xfprintf(out,"%d\n",elm->id);

        xfclose(out);
     */

    if (OptGetBool("Transport", "Transport_on", "no") == true) {
    	problem->otransport->convection();
    }

    /*
        if (OptGetBool("Transport",  "Reactions", "no") == true) {
            read_reaction_list(transport);
        }
*/
        //if (rank == 0) {

        //}

        // TODO: there is an uncoditioned jump in open_temp_files
        // also this function should be moved to btc.*
        // btc should be documented and have an clearly defined interface
        // not strictly dependent on Transport.
        //btc_check(transport);



        /*
                if(problem->cross_section == true)
                {
                    elect_cross_section_element(problem);
                    output_transport_init_CS(problem);
                    output_transport_time_CS(problem, 0 * problem->time_step);
                }
         */
}
//-----------------------------------------------------------------------------
// vim: set cindent:
//-----------------------------------------------------------------------------
#endif
#if 0

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR UNSTEADY SATURATED FLOW
 */
void main_compute_mh_density(struct Problem *problem)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    int i, j, dens_step, n_step, frame = 0, rank;
    double save_step, stop_time; // update_dens_time
    char statuslog[255];
    FILE *log;
    OutputTime *output_time = NULL;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

/*
    transport_output_init(problem->otransport->transport_out_fname);
    transport_output(problem->otransport->out_conc,problem->otransport->substance_name ,problem->otransport->n_substances, 0.0, ++frame,problem->otransport->transport_out_fname);
*/

    if(rank == 0) {
        output_time = new OutputTime(mesh, problem->otransport->transport_out_fname);
    }

    //save_step = problem->save_step;
    //stop_time = problem->stop_time;
    //trans->update_dens_time = problem->save_step / (ceil(problem->save_step / trans->dens_step));
    //dens_step = (int) ceil(problem->stop_time / trans->update_dens_time);
    //n_step = (int) (problem->save_step / trans->update_dens_time);

    // DF problem - I don't understend to this construction !!!
    //problem->save_step = problem->stop_time = trans->update_dens_time;


    /*



    //------------------------------------------------------------------------------
    //      Status LOG head
    //------------------------------------------------------------------------------
    //sprintf( statuslog,"%s.txt",problem->log_fname);
    sprintf(statuslog, "density_log.txt");
    log = xfopen(statuslog, "wt");

    xfprintf(log, "Stop time = %f (%f) \n", trans->update_dens_time * dens_step, stop_time);
    xfprintf(log, "Save step = %f \n", save_step);
    xfprintf(log, "Density  step = %f (%d) \n\n", trans->update_dens_time, trans->dens_step);
    xfprintf(log, "Time\t Iteration number\n");
    //------------------------------------------------------------------------------

    for (i = 0; i < dens_step; i++) {
        xprintf(Msg, "dens step %d \n", i);
        save_time_step_C(problem);

        for (j = 0; j < trans->max_dens_it; j++) {
            xprintf(Msg, "dens iter %d \n", j);
            save_restart_iteration_H(problem);
            //restart_iteration(problem);
            //calculation_mh(problem);
            //problem->water=new DarcyFlowMH(*mesh);
            //problem->water->solve();
            restart_iteration_C(problem);
            //postprocess(problem);
            convection(trans, output_time);

            if (trans->dens_implicit == 0) {
                xprintf(Msg, "no density iterations (explicit)", j);
                break;
            }
            if (compare_dens_iter(problem) && (j > 0)) {
                break; //at least one repeat of iteration is necessary to update both conc and pressure
            }
        }

        xprintf(Msg, "step %d finished at %d density iterations\n", i, j);
        xfprintf(log, "%f \t %d\n", (i + 1) * trans->update_dens_time, j); // Status LOG
    }

    if(rank == 0) {
        delete output_time;
    }

    xfclose(log); */
//}
#endif

