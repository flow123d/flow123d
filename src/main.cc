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
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include <iostream>
#include <fstream>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>

#include "main.h"

#include "rev_num.h"

/// named version of the program
#define _PROGRAM_VERSION_   "0.0.0"

#ifndef _PROGRAM_REVISION_
    #define _PROGRAM_REVISION_ "(unknown revision)"
#endif

#ifndef _PROGRAM_BRANCH_
    #define _PROGRAM_BRANCH_ "(unknown branch)"
#endif

#ifndef FLOW123D_COMPILER_FLAGS_
    #define FLOW123D_COMPILER_FLAGS_ "(unknown compiler flags)"
#endif


namespace it = Input::Type;

// this should be part of a system class containing all support information
it::Record & Application::get_input_type() {
    static it::Record type = it::Record("Root", "Root record of JSON input for Flow123d.")
    .declare_key("problem", CouplingBase::get_input_type(), it::Default::obligatory(),
    		"Simulation problem to be solved.")
    .declare_key("pause_after_run", it::Bool(), it::Default("false"),
    		"If true, the program will wait for key press before it terminates.")
	.close();

    return type;
}



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
#ifdef FLOW123D_HAVE_PYTHON
    PythonLoader::initialize(argv[0]);
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


Input::Type::RevNumData Application::get_rev_num_data() {
	Input::Type::RevNumData rev_num_data;
	rev_num_data.version = string(_VERSION_NAME_);
	rev_num_data.revision = string(_GIT_REVISION_);
	rev_num_data.branch = string(_GIT_BRANCH_);
	rev_num_data.url = string(_GIT_URL_);

	return rev_num_data;
}


void Application::display_version() {
    // Say Hello
    // make strings from macros in order to check type
	Input::Type::RevNumData rev_num_data = this->get_rev_num_data();
    string build = string(__DATE__) + ", " + string(__TIME__) + " flags: " + string(FLOW123D_COMPILER_FLAGS_);


    xprintf(Msg, "This is Flow123d, version %s revision: %s\n", rev_num_data.version.c_str(), rev_num_data.revision.c_str());
    xprintf(Msg,
    	 "Branch: %s\n"
		 "Build: %s\n"
		 "Fetch URL: %s\n",
		 rev_num_data.branch.c_str(), build.c_str() , rev_num_data.url.c_str() );
    Profiler::instance()->set_program_info("Flow123d", rev_num_data.version, rev_num_data.branch, rev_num_data.revision, build);
}



Input::Record Application::read_input() {
   if (main_input_filename_ == "") {
        cout << "Usage error: The main input file has to be specified through -s parameter.\n\n";
        cout << program_arguments_desc_ << "\n";
        exit( exit_failure );
    }

    // read main input file
    FilePath fpath(main_input_filename_, FilePath::FileType::input_file);
    try {
    	Input::ReaderToStorage json_reader(fpath, get_input_type() );
        root_record = json_reader.get_root_interface<Input::Record>();
    } catch (Input::ReaderToStorage::ExcInputError &e ) {
      e << Input::ReaderToStorage::EI_File(fpath); throw;
    } catch (Input::ReaderToStorage::ExcNotJSONFormat &e) {
      e << Input::ReaderToStorage::EI_File(fpath); throw;
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
        ("JSON_machine", po::value< string >(), "Writes full structure of the main input file as a valid CON file into given file")
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
    /*if (vm.count("full_doc")) {
        Input::Type::TypeBase::lazy_finish();
        Input::Type::OutputText type_output(&get_input_type());
        type_output.set_filter(":Field:.*");
        cout << type_output;
        exit( exit_output );
    }*/

    // if there is "JSON_machine" option
    if (vm.count("JSON_machine")) {
        // write ist to json file
        string json_filename = vm["JSON_machine"].as<string>();
        ofstream json_stream(json_filename);
        // check open operation
        if (json_stream.fail()) {
    		cerr << "Failed to open file '" << json_filename << "'" << endl;
        } else {
	        Input::Type::TypeBase::lazy_finish();
	        json_stream << Input::Type::OutputJSONMachine( this->get_rev_num_data() );
	        json_stream.close();
        }
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

	display_version();

    Input::Record i_rec = read_input();


    {
        using namespace Input;

        // get main input record handle

        // should flow123d wait for pressing "Enter", when simulation is completed
        sys_info.pause_after_run = i_rec.val<bool>("pause_after_run");
        // read record with problem configuration
        Input::AbstractRecord i_problem = i_rec.val<AbstractRecord>("problem");

        if (i_problem.type() == HC_ExplicitSequential::get_input_type() ) {

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
        // log profiler data to this stream
        Profiler::instance()->output (PETSC_COMM_WORLD);

        // call python script which transforms json file at given location
        // Profiler::instance()->transform_profiler_data (".csv", "CSVFormatter");
        Profiler::instance()->transform_profiler_data (".txt", "SimpleTableFormatter");

        // finally uninitialize
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
