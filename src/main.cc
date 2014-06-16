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

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/python_loader.hh"
#include "coupling/hc_explicit_sequential.hh"
#include "input/input_type.hh"
#include "input/type_output.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"
//#include "io/output.h"

#include <iostream>
#include <fstream>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem.hpp>

#include "main.h"
//#include "io/read_ini.h"

#include "rev_num.h"

/// named version of the program
#define _PROGRAM_VERSION_   "0.0.0"

#ifndef _PROGRAM_REVISION_
    #define _PROGRAM_REVISION_ "(unknown revision)"
#endif

#ifndef _PROGRAM_BRANCH_
    #define _PROGRAM_BRANCH_ "(unknown branch)"
#endif

#ifndef _COMPILER_FLAGS_
    #define _COMPILER_FLAGS_ "(unknown compiler flags)"
#endif

//static void main_convert_to_output();


namespace it = Input::Type;

// this should be part of a system class containing all support information
//static Record system_rec("System", "Record with general support data.");
//system_rec.finish();
it::Record Application::input_type
    = it::Record("Root", "Root record of JSON input for Flow123d.")
    //main_rec.declare_key("system", system_rec, "");
    .declare_key("problem", CouplingBase::input_type, it::Default::obligatory(),
    		"Simulation problem to be solved.")
    .declare_key("pause_after_run", it::Bool(), it::Default("false"),
    		"If true, the program will wait for key press before it terminates.");
//    .declare_key("output_streams", it::Array( OutputTime::input_type ),
//    		"Array of formated output streams to open.");



Application::Application( int argc,  char ** argv)
: ApplicationBase(argc, argv),
  main_input_dir_("."),
  main_input_filename_(""),
  passed_argc_(0),
  passed_argv_(0),
  use_profiler(true)
{
    // initialize python stuff if we have
    // nonstandard python home (release builds)
    std::cout << "Application constructor" << std::endl;
#ifdef HAVE_PYTHON
#ifdef PYTHON_HOME
    PythonLoader::initialize(argv[0]);
#endif
#endif

}


void Application::split_path(const string& path, string& directory, string& file_name) {

    size_t delim_pos=path.find_last_of(DIR_DELIMITER);
    if (delim_pos < string::npos) {

        // It seems, that there is some path in fname ... separate it
        directory =path.substr(0,delim_pos);
        file_name =path.substr(delim_pos+1); // till the end
    } else {
        directory = ".";
        file_name = path;
    }
}

void Application::display_version() {
    // Say Hello
    // make strings from macros in order to check type
    string version(_VERSION_NAME_);
    string revision(_GIT_REVISION_);
    string branch(_GIT_BRANCH_);
    string url(_GIT_URL_);
    string build = string(__DATE__) + ", " + string(__TIME__) + " flags: " + string(_COMPILER_FLAGS_);
    

    xprintf(Msg, "This is Flow123d, version %s revision: %s\n", version.c_str(), revision.c_str());
    xprintf(Msg,
    	 "Branch: %s\n"
		 "Build: %s\n"
		 "Fetch URL: %s\n",
		 branch.c_str(), build.c_str() , url.c_str() );
    Profiler::instance()->set_program_info("Flow123d", version, branch, revision, build);
}



Input::Record Application::read_input() {
   if (main_input_filename_ == "") {
        cout << "Usage error: The main input file has to be specified through -s parameter.\n\n";
        cout << program_arguments_desc_ << "\n";
        exit( exit_failure );
    }
    
    // read main input file
    string fname = main_input_dir_ + DIR_DELIMITER + main_input_filename_;
    std::ifstream in_stream(fname.c_str());
    if (! in_stream) {
        xprintf(UsrErr, "Can not open main input file: '%s'.\n", fname.c_str());
    }
    try {
    	Input::JSONToStorage json_reader(in_stream, input_type );
        root_record = json_reader.get_root_interface<Input::Record>();
    } catch (Input::JSONToStorage::ExcInputError &e ) {
      e << Input::JSONToStorage::EI_File(fname); throw;
    } catch (Input::JSONToStorage::ExcNotJSONFormat &e) {
      e << Input::JSONToStorage::EI_File(fname); throw;
    }  
    
    return root_record;
}




void Application::parse_cmd_line(const int argc, char ** argv) {
	namespace po = boost::program_options;


    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("solve,s", po::value< string >(), "Main input file to solve.")
        ("input_dir,i", po::value< string >()->default_value("input"), "Directory for the ${INPUT} placeholder in the main input file.")
        ("output_dir,o", po::value< string >()->default_value("output"), "Directory for all produced output files.")
        ("log,l", po::value< string >()->default_value("flow123"), "Set base name for log files.")
        ("version", "Display version and build information and exit.")
        ("no_log", "Turn off logging.")
        ("no_profiler", "Turn off profiler output.")
        ("full_doc", "Prints full structure of the main input file.")
        ("JSON_template", "Prints description of the main input file as a valid CON file.")
        ("latex_doc", "Prints description of the main input file in Latex format using particular macros.")
    	("JSON_machine", "Prints full structure of the main input file as a valid CON file.")
        ("petsc_redirect", po::value<string>(), "Redirect all PETSc stdout and stderr to given file.");

    ;

    // parse the command line
    po::variables_map vm;
    po::parsed_options parsed = po::basic_command_line_parser<char>(argc, argv).options(desc).allow_unregistered().run();
    po::store(parsed, vm);
    po::notify(vm);

    // get unknown options
    vector<string> to_pass_further = po::collect_unrecognized(parsed.options, po::include_positional);
    passed_argc_ = to_pass_further.size();
    passed_argv_ = new char * [passed_argc_+1];

    // first copy the program executable in argv[0]
    int arg_i=0;
    if (argc > 0) passed_argv_[arg_i++] = xstrcpy( argv[0] );

    for(int i=0; i < passed_argc_; i++) {
        passed_argv_[arg_i++] = xstrcpy( to_pass_further[i].c_str() );
    }
    passed_argc_ = arg_i;

    // if there is "help" option
    if (vm.count("help")) {
        cout << desc << "\n";
        exit( exit_output );
    }

    if (vm.count("version")) {
    	display_version();
    	exit( exit_output );
    }

    // if there is "full_doc" option
    if (vm.count("full_doc")) {
        Input::Type::TypeBase::lazy_finish();
        Input::Type::OutputText type_output(&input_type);
        type_output.set_filter(":Field:.*");
        cout << type_output;
        exit( exit_output );
    }

    if (vm.count("JSON_template")) {
        Input::Type::TypeBase::lazy_finish();
        cout << Input::Type::OutputJSONTemplate(&input_type);
        exit( exit_output );
    }

    if (vm.count("latex_doc")) {
        Input::Type::TypeBase::lazy_finish();
        Input::Type::OutputLatex type_output(&input_type);
        type_output.set_filter("");
        cout << type_output;
        exit( exit_output );
    }

    if (vm.count("JSON_machine")) {
        Input::Type::TypeBase::lazy_finish();
        cout << Input::Type::OutputJSONMachine(&input_type);
        exit( exit_output );
    }

    if (vm.count("petsc_redirect")) {
        this->petsc_redirect_file_ = vm["petsc_redirect"].as<string>();
    }

    // if there is "solve" option
    if (vm.count("solve")) {
        string input_filename = vm["solve"].as<string>();
        split_path(input_filename, main_input_dir_, main_input_filename_);
    } 

    // possibly turn off profilling
    if (vm.count("no_profiler")) use_profiler=false;

    string input_dir;
    string output_dir;
    if (vm.count("input_dir")) {
        input_dir = vm["input_dir"].as<string>();
    }
    if (vm.count("output_dir")) {
            output_dir = vm["output_dir"].as<string>();
    }

    // assumes working directory "."
    FilePath::set_io_dirs(".", main_input_dir_, input_dir, output_dir );

    if (!boost::filesystem::is_directory(output_dir)) {
    	boost::filesystem::create_directory(output_dir);
    }

    if (vm.count("log")) {
        this->log_filename_ = vm["log"].as<string>();
    }

    if (vm.count("no_log")) {
        this->log_filename_="//";     // override; do not open log files
    }

    ostringstream tmp_stream(program_arguments_desc_);
    tmp_stream << desc;
    // TODO: catch specific exceptions and output usage messages
}





void Application::run() {
    //use_profiler=true;


    display_version();

    Input::Record i_rec = read_input();


    {
        using namespace Input;

        // get main input record handle

        // should flow123d wait for pressing "Enter", when simulation is completed
        sys_info.pause_after_run = i_rec.val<bool>("pause_after_run");
        // read record with problem configuration
        Input::AbstractRecord i_problem = i_rec.val<AbstractRecord>("problem");

        if (i_problem.type() == HC_ExplicitSequential::input_type ) {

            HC_ExplicitSequential *problem = new HC_ExplicitSequential(i_problem);

            // run simulation
            problem->run_simulation();

            delete problem;
        } else {
            xprintf(UsrErr,"Problem type not implemented.");
        }

    }
}




void Application::after_run() {
	if (sys_info.pause_after_run) {
        printf("\nPress <ENTER> for closing the window\n");
        getchar();
    }
}




Application::~Application() {
    if (use_profiler && Profiler::is_initialized()) {
        Profiler::instance()->output(PETSC_COMM_WORLD);
        Profiler::uninitialize();
    }
}


//=============================================================================

/**
 *  FUNCTION "MAIN"
 */
int main(int argc, char **argv) {
    try {
        Application app(argc, argv);
        app.init(argc, argv);
    } catch (std::exception & e) {
        std::cerr << e.what();
        return ApplicationBase::exit_failure;
    } catch (...) {
        std::cerr << "Unknown exception" << endl;
        return ApplicationBase::exit_failure;
    }

    // Say Goodbye
    return ApplicationBase::exit_success;
}






























/**
 * FUNCTION "MAIN" FOR CONVERTING FILES TO POS
 */
/*void main_convert_to_output() {
    // TODO: implement output of input data fields
    // Fields to output:
    // 1) volume data (simple)
    //    sources (Darcy flow and transport), initial condition, material id, partition id
    // 2) boundary data (needs "virtual fractures" in output mesh)
    //    flow and transport bcd

    xprintf(Err, "Not implemented yet in this version\n");
}*/
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

