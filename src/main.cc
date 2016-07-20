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
 * @file    main.cc
 * @brief   This file should contain only creation of Application object.
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
#include <boost/regex.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem.hpp>

#include "main.h"

#include "rev_num.h"

/// named version of the program
//#define _PROGRAM_VERSION_   "0.0.0"

//#ifndef _PROGRAM_REVISION_
//    #define _PROGRAM_REVISION_ "(unknown revision)"
//#endif

//#ifndef _PROGRAM_BRANCH_
//    #define _PROGRAM_BRANCH_ "(unknown branch)"
//#endif

#ifndef FLOW123D_COMPILER_FLAGS_
    #define FLOW123D_COMPILER_FLAGS_ "(unknown compiler flags)"
#endif


namespace it = Input::Type;

// this should be part of a system class containing all support information
it::Record & Application::get_input_type() {
    static it::Record type = it::Record("Root", "Root record of JSON input for Flow123d.")
    .declare_key("flow123d_version", it::String(), it::Default::obligatory(),
            "Version of Flow123d for which the input file was created."
            "Flow123d only warn about version incompatibility. "
            "However, external tools may use this information to provide conversion "
            "of the input file to the structure required by another version of Flow123d.")
    .declare_key("problem", CouplingBase::get_input_type(), it::Default::obligatory(),
    		"Simulation problem to be solved.")
    .declare_key("pause_after_run", it::Bool(), it::Default("false"),
    		"If true, the program will wait for key press before it terminates.")
	.close();

    return type;
}



Application::Application( int argc,  char ** argv)
: ApplicationBase(argc, argv),
  main_input_filename_(""),
  passed_argc_(0),
  passed_argv_(0),
  use_profiler(true),
  yaml_balance_output_(false)
{
    // initialize python stuff if we have
    // nonstandard python home (release builds)
#ifdef FLOW123D_HAVE_PYTHON
    PythonLoader::initialize(argv[0]);
#endif

}


Input::Type::RevNumData Application::get_rev_num_data() {
	Input::Type::RevNumData rev_num_data;

	rev_num_data.version = string(FLOW123D_VERSION_NAME_);
	rev_num_data.revision = string(FLOW123D_GIT_REVISION_);
	rev_num_data.branch = string(FLOW123D_GIT_BRANCH_);
	rev_num_data.url = string(FLOW123D_GIT_URL_);

	return rev_num_data;
}


void Application::display_version() {
    // Say Hello
    // make strings from macros in order to check type
	Input::Type::RevNumData rev_num_data = this->get_rev_num_data();
    string build = string(__DATE__) + ", " + string(__TIME__)
            + " flags: " + string(FLOW123D_COMPILER_FLAGS_);


    xprintf(Msg, "This is Flow123d, version %s revision: %s\n",
            rev_num_data.version.c_str(),
            rev_num_data.revision.c_str());
    xprintf(Msg,
    	 "Branch: %s\n"
		 "Build: %s\n"
		 "Fetch URL: %s\n",
		 rev_num_data.branch.c_str(), build.c_str() , rev_num_data.url.c_str() );
    Profiler::instance()->set_program_info("Flow123d",
            rev_num_data.version, rev_num_data.branch, rev_num_data.revision, build);
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
		("petsc_redirect", po::value<string>(), "Redirect all PETSc stdout and stderr to given file.")
		("yaml_balance", "Redirect balance output to YAML format too (simultaneously with the selected balance output format).");

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
            // create the root Record
            it::Record root_type = get_input_type();
            Input::Type::TypeBase::lazy_finish();
            json_stream << Input::Type::OutputJSONMachine( root_type, this->get_rev_num_data() );
            json_stream.close();
        }
        exit( exit_output );
    }

    if (vm.count("petsc_redirect")) {
        this->petsc_redirect_file_ = vm["petsc_redirect"].as<string>();
    }

    // if there is "solve" option
    string input_filename = ".";
    if (vm.count("solve")) {
        input_filename = vm["solve"].as<string>();
    }

    // possibly turn off profilling
    if (vm.count("no_profiler")) use_profiler=false;

    // preserves output of balance in YAML format
    if (vm.count("yaml_balance")) yaml_balance_output_=true;

    string input_dir;
    string output_dir;
    if (vm.count("input_dir")) {
        input_dir = vm["input_dir"].as<string>();
    }
    if (vm.count("output_dir")) {
            output_dir = vm["output_dir"].as<string>();
    }

    // assumes working directory "."
    main_input_filename_ = FilePath::set_dirs_from_input(input_filename, input_dir, output_dir );

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
    START_TIMER("Application::run");
	display_version();

    START_TIMER("Read Input");
    // get main input record handle
    Input::Record i_rec = read_input();
    END_TIMER("Read Input");

    {
        using namespace Input;
        // check input file version against the version of executable
        boost::regex version_re("([^.]*)[.]([^.]*)[.]([^.]*)");
        boost::smatch match;
        std::string version(FLOW123D_VERSION_NAME_);
        vector<string> ver_fields(3);
        if ( boost::regex_match(version, match, version_re) ) {
            ver_fields[0]=match[1];
            ver_fields[1]=match[2];
            ver_fields[2]=match[3];
        } else {
        	OLD_ASSERT(1, "Bad Flow123d version format: %s\n", version.c_str() );
        }

        std::string input_version = i_rec.val<string>("flow123d_version");
        vector<string> iver_fields(3);
        if ( boost::regex_match(input_version, match, version_re) ) {
            iver_fields[0]=match[1];
            iver_fields[1]=match[2];
            iver_fields[2]=match[3];
        } else {
            THROW( ExcVersionFormat() << EI_InputVersionStr(input_version) );
        }

        if ( iver_fields[0] != ver_fields[0] || iver_fields[1] > ver_fields[1] ) {
            xprintf(Warn, "Input file with version: '%s' is no compatible with the program version: '%s' \n",
                    input_version.c_str(), version.c_str());
        }

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
	// remove balance output files in YAML format if "yaml_balance" option is not set
	if ( (sys_info.my_proc==0) && !yaml_balance_output_ ) {
		boost::filesystem::path mass_file( string(FilePath("mass_balance.yaml", FilePath::output_file)) );
		boost::filesystem::path water_file( string(FilePath("water_balance.yaml", FilePath::output_file)) );
		boost::filesystem::path energy_file( string(FilePath("energy_balance.yaml", FilePath::output_file)) );

		if (boost::filesystem::exists(mass_file)) {
		    boost::filesystem::remove(mass_file);
		}
		if (boost::filesystem::exists(water_file)) {
		    boost::filesystem::remove(water_file);
		}
		if (boost::filesystem::exists(energy_file)) {
		    boost::filesystem::remove(energy_file);
		}
	}

    if (use_profiler) {
        if (petsc_initialized) {
            // log profiler data to this stream
            Profiler::instance()->output (PETSC_COMM_WORLD);
        } else {
            Profiler::instance()->output();
        }

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
        _LOG( Logger::MsgType::error ) << e.what();
        return ApplicationBase::exit_failure;
    } catch (...) {
        _LOG( Logger::MsgType::error ) << "Unknown exception" << endl;
        return ApplicationBase::exit_failure;
    }

    // Say Goodbye
    return ApplicationBase::exit_success;
}
