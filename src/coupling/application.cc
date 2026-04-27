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
 * @file    application_base.cc
 * @brief   
 */

#include "application.hh"

#include "system/sys_profiler.hh"
#include "system/logger_options.hh"
#include "system/file_path.hh"
#include "system/system.hh"
#include <signal.h>
#include <iostream>

#ifdef FLOW123D_HAVE_PETSC
//#include <petsc.h>
#include <petscsys.h>
#include <petsc/private/petscimpl.h> /* to gain access to the private PetscVFPrintf */
#endif

#ifdef FLOW123D_HAVE_PERMON
#include <permonsys.h>
#endif

#include <string.h>                                    // for strsignal

#include <iostream>                                    // for cout
#include <sstream>                                     // for operator<<, endl
#include "mpi.h"                                       // for MPI_Comm_size
#include "petscerror.h"                                // for CHKERRQ, Petsc...
#include "system/exc_common.hh"                        // for ExcAssertMsg
#include "system/asserts.hh"                           // for ASSERT_PERMANENT, msg
#include "system/logger.hh"                            // for Logger, operat...
#include "system/system.hh"                            // for SystemInfo




/// Function that catches all program signals.
/// Note: context variable required by PETSc function PetscPushSignalHandler
PetscErrorCode petsc_signal_handler(int signal, FMT_UNUSED void *context)
{
  if (signal == SIGINT) {
      cout << "SIGINT\n";
  }
  if (signal == SIGFPE ||   // FPE: Floating Point Exception,probably divide by zero
      signal == SIGILL ||   // Illegal instruction: Likely due to memory corruption
      signal == SIGPIPE ||  // Broken Pipe: Likely while reading or writing to a socket
      signal == SIGSEGV )   // SEGV: Segmentation Violation, probably memory access out of range
  {
      // Signals handled by us.
      THROW( ExcSignal() << EI_Signal(signal) << EI_SignalName(strsignal(signal)) );
  } else {
      return PetscSignalHandlerDefault(signal,(void*)0);
  }
  return 0;
}

void system_signal_handler(int signal) {
    petsc_signal_handler(signal, nullptr);
}


Application::Application()
: log_filename_(""),
  signal_handler_off_(false),
  problem_(nullptr),
  main_input_filename_(""),
  //passed_argc_(0),
  //passed_argv_(0),
  use_profiler(true),
  memory_monitoring(false),
  profiler_path(""),
  yaml_balance_output_(false)

{
    // initialize python stuff if we have
    // nonstandard python home (release builds)
    PythonLoader::initialize();

}




bool Application::petsc_initialized = false;
bool Application::permon_initialized = false;


void Application::system_init( MPI_Comm comm, const string &log_filename ) {
    int ierr;

    sys_info.comm=comm;


    //Xio::init(); //Initialize XIO library

    // TODO : otevrit docasne log file jeste pred ctenim vstupu (kvuli zachyceni chyb), po nacteni dokoncit
    // inicializaci systemu

    ierr=MPI_Comm_rank(comm, &(sys_info.my_proc));
    ierr+=MPI_Comm_size(comm, &(sys_info.n_proc));
    LoggerOptions::get_instance().setup_mpi(comm);
    ASSERT_PERMANENT( ierr == MPI_SUCCESS ).error("MPI not initialized.\n");

    // determine logfile name or switch it off
    stringstream log_name;

    if ( log_filename == "//" ) {
    	// -l option without given name -> turn logging off
    	sys_info.log=NULL;
    	LoggerOptions::get_instance().set_log_file("");
    } else	{
    	// construct full log name
    	//log_name << log_filename <<  "." << sys_info.my_proc << ".old.log";

    	//sys_info.log_fname = FilePath(log_name.str(), FilePath::output_file);
    	//sys_info.log=xfopen(sys_info.log_fname.c_str(),"wt");

    	LoggerOptions::get_instance().set_log_file(log_filename);
    }

    sys_info.verbosity=0;
    sys_info.pause_after_run=0;
}


FILE *Application::petsc_output_ =NULL;

#ifdef FLOW123D_HAVE_PETSC
PetscErrorCode Application::petscvfprintf(FILE *fd, const char format[], va_list Argp) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (fd != stdout && fd != stderr) { /* handle regular files */
    ierr = PetscVFPrintfDefault(fd,format,Argp); CHKERRQ(ierr);
  } else {
    const int buf_size = 65000;
    char buff[65000];
    size_t length;
    ierr = PetscVSNPrintf(buff,buf_size,format,&length,Argp);CHKERRQ(ierr);

    /* now send buff to whatever stream or whatever you want */
    fwrite(buff, sizeof(char), length, petsc_output_);
  }
  PetscFunctionReturn(0);
}
#endif


void Application::petsc_initialize(int argc, char ** argv) {
#ifdef FLOW123D_HAVE_PETSC
    if (petsc_redirect_file_ != "") {
        petsc_output_ = fopen(petsc_redirect_file_.c_str(), "w");
        if (! petsc_output_)
            THROW(FilePath::ExcFileOpen() << FilePath::EI_Path(petsc_redirect_file_));
        PetscVFPrintf = this->petscvfprintf;
    }


    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    if (! signal_handler_off_) {
        // PETSc do not catch SIGINT, but someone on the way does, we try to fix it.
        signal(SIGINT, system_signal_handler);
        PetscPushSignalHandler(petsc_signal_handler, nullptr);
    }

    int mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MessageOut() << "MPI size: " << mpi_size << std::endl;
#endif
}



int Application::petcs_finalize() {
#ifdef FLOW123D_HAVE_PETSC
	if ( petsc_initialized )
	{
		PetscErrorCode ierr=0;

		ierr = PetscFinalize(); CHKERRQ(ierr);

		if (petsc_output_) fclose(petsc_output_);

		petsc_initialized = false;

		return ierr;
	}
#endif

	return 0;
}


void Application::permon_initialize(int argc, char ** argv) {
#ifdef FLOW123D_HAVE_PERMON
    PermonInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
#endif
}

int Application::permon_finalize() {
#ifdef FLOW123D_HAVE_PERMON
	if ( permon_initialized )
	{
		PetscErrorCode ierr=0;

		ierr = PermonFinalize(); CHKERRQ(ierr);

		permon_initialized = false;

		return ierr;
	}
#endif

	return 0;
}


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

Input::Type::RevNumData get_rev_num_data() {
	Input::Type::RevNumData rev_num_data;

	rev_num_data.version = string(FLOW123D_VERSION_NAME_);
	rev_num_data.revision = string(FLOW123D_GIT_REVISION_);
	rev_num_data.branch = string(FLOW123D_GIT_BRANCH_);
	rev_num_data.url = string(FLOW123D_GIT_URL_);

	return rev_num_data;
}


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







void Application::display_version() {
    // Say Hello
    // make strings from macros in order to check type
	Input::Type::RevNumData rev_num_data = get_rev_num_data();


	string build = string(__DATE__) + ", " + string(__TIME__)
            + " flags: " + string(FLOW123D_COMPILER_FLAGS_);


    MessageOut().fmt("This is Flow123d, version {} commit: {}\n",
            rev_num_data.version, rev_num_data.revision);
    MessageOut().fmt("Branch: {}\nBuild: {}\nFetch URL: {}\n",
		 rev_num_data.branch, build, rev_num_data.url );

}



Input::Record Application::read_input() {
   if (main_input_filename_ == "") {
        cout << "Usage error: The main input file has to be specified through -s parameter.\n\n";
        cout << program_arguments_desc_ << "\n";
        exit( exit_failure );
    }

    // read main input file
    FilePath fpath(main_input_filename_, FilePath::FileType::input_file);

    Input::ReaderToStorage json_reader(fpath, get_input_type() );
    root_record = json_reader.get_root_interface<Input::Record>();

    return root_record;
}



void Application::parse_cmd_line(const int argc, char ** argv) {
	namespace po = boost::program_options;



    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("solve,s", po::value< string >(), "Main input file to solve.")
        ("input_dir,i", po::value< string >()->default_value("input"), "Directory for the $INPUT_DIR$ placeholder in the main input file.")
        ("output_dir,o", po::value< string >()->default_value("output"), "Directory for all produced output files.")
        ("log,l", po::value< string >()->default_value("flow123"), "Set base name for log files.")
        ("version", "Display version and build information and exit.")
        ("no_log", "Turn off logging.")
        ("no_signal_handler", "Turn off signal handling. Useful for debugging with valgrind.")
        ("no_profiler,no-profiler", "Turn off profiler output.")
        ("profiler_path,profiler-path", po::value< string >(), "Path to the profiler file")
        ("memory_monitoring,memory-monitoring", "Switch on memory monitoring.")
        ("input_format", po::value< string >(), "Writes full structure of the main input file into given file.")
		("petsc_redirect", po::value<string>(), "Redirect all PETSc stdout and stderr to given file.")
		("yaml_balance", "Redirect balance output to YAML format too (simultaneously with the selected balance output format).");



    // Can not use positional arguments together with PETSC options.
    // Use our own solution trying to use the first unrecognized option as the main input file.

    // parse the command line
    po::variables_map vm;
    auto parser = po::basic_command_line_parser<char>(argc, argv)
            .options(desc)
            .allow_unregistered();
    po::parsed_options parsed = parser.run();
    po::store(parsed, vm);
    po::notify(vm);


    // get unknown options
    vector<string> to_pass_further = po::collect_unrecognized(parsed.options, po::include_positional);



    /*
    passed_argc_ = to_pass_further.size();
    passed_argv_ = new char * [passed_argc_+1];

    // first copy the program executable in argv[0]
    int arg_i=0;
    if (argc > 0) passed_argv_[arg_i++] = xstrcpy( argv[0] );

    for(int i=0; i < passed_argc_; i++) {
        passed_argv_[arg_i++] = xstrcpy( to_pass_further[i].c_str() );
    }
    passed_argc_ = arg_i;
    */

    // possibly turn off profilling
    if (vm.count("no_profiler")) {
        use_profiler=false;
    }

    if (vm.count("profiler_path")) {
        profiler_path = vm["profiler_path"].as<string>();
    }

    if (vm.count("memory_monitoring")) {
        memory_monitoring=true;
    }

    // if there is "help" option
    if (vm.count("help")) {
        display_version();
        cout << endl;
        cout << "Usage:" << endl;
        cout << "   flow123d -s <main_input>.yaml <other options> <PETSC options>" << endl;
        cout << "   flow123d <main_input>.yaml <other options> <PETSC options>" << endl;
        cout << desc << "\n";
        THROW(ExcNoRunOption());
    }


    if (vm.count("version")) {
    	display_version();
    	THROW(ExcNoRunOption());
    }



    // if there is "input_format" option
    if (vm.count("input_format")) {
        // write ist to json file
        ofstream json_stream;
        FilePath(vm["input_format"].as<string>(), FilePath::output_file).open_stream(json_stream);
        // create the root Record
        it::Record root_type = get_input_type();
        root_type.finish();
        Input::Type::TypeBase::delete_unfinished_types();
        json_stream << Input::Type::OutputJSONMachine( root_type, get_rev_num_data() );
        json_stream.close();
        THROW(ExcNoRunOption());
    }


    if (vm.count("petsc_redirect")) {
        this->petsc_redirect_file_ = vm["petsc_redirect"].as<string>();
    }

    if (vm.count("no_signal_handler")) {
        this->signal_handler_off_ = true;
    }

    // if there is "solve" option
    string input_filename = "";

    // check for positional main input file
    if (to_pass_further.size()) {
        string file_candidate = to_pass_further[0];
        if (file_candidate[0] != '-') {
            // pop the first option
            input_filename = file_candidate;
            to_pass_further.erase(to_pass_further.begin());
        }
    }


    if (vm.count("solve")) {
        input_filename = vm["solve"].as<string>();
    }


    if (input_filename == "")
        THROW(ExcMessage() << EI_Message("Main input file not specified (option -s)."));

    // preserves output of balance in YAML format
    if (vm.count("yaml_balance")) Balance::set_yaml_output();

    string input_dir;
    string output_dir;
    if (vm.count("input_dir")) {
        input_dir = vm["input_dir"].as<string>();
    }
    if (vm.count("output_dir")) {
            output_dir = vm["output_dir"].as<string>();
    }



    // assumes working directory "."
    try {
        main_input_filename_ = FilePath::set_dirs_from_input(input_filename, input_dir, output_dir );
    } catch (FilePath::ExcMkdirFail &e) {
        use_profiler = false; // avoid profiler output
        throw e;
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

/**
 * Contains basic structure of application (initialization, run and finalization).
 * Method is call after constructor and allows to call virtual methods.
 */

void Application::init(int argc, char ** argv) {
    // parse our own command line arguments, leave others for PETSc

	this->parse_cmd_line(argc, argv);

	string build = string(__DATE__) + ", " + string(__TIME__)
            + " flags: " + string(FLOW123D_COMPILER_FLAGS_);

	Input::Type::RevNumData rev_num_data = get_rev_num_data();
    Profiler::instance()->set_program_info("Flow123d",
            rev_num_data.version, rev_num_data.branch, rev_num_data.revision, build);

    if (use_profiler & memory_monitoring)
        Profiler::set_memory_monitoring(memory_monitoring);

    armadillo_setup(); // set catching armadillo exceptions and reporting stacktrace

	this->petsc_initialize(argc, argv);
	petsc_initialized = true;

	this->permon_initialize(argc, argv);
	permon_initialized = true;

    this->system_init(PETSC_COMM_WORLD, log_filename_); // Petsc, open log, read ini file


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
        std::regex version_re("([^.]*)[.]([^.]*)[.]([^.]*)");
        std::smatch match;
        std::string version(FLOW123D_VERSION_NAME_);
        vector<string> ver_fields(3);
        if ( std::regex_match(version, match, version_re) ) {
            ver_fields[0]=match[1];
            ver_fields[1]=match[2];
            ver_fields[2]=match[3];
        } else {
        	ASSERT_PERMANENT(0)(version).error("Bad Flow123d version format\n");
        }

        std::string input_version = i_rec.val<string>("flow123d_version");
        vector<string> iver_fields(3);
        if ( std::regex_match(input_version, match, version_re) ) {
            iver_fields[0]=match[1];
            iver_fields[1]=match[2];
            iver_fields[2]=match[3];
        } else {
            THROW( ExcVersionFormat() << EI_InputVersionStr(input_version) );
        }

        if ( iver_fields[0] != ver_fields[0] || iver_fields[1] > ver_fields[1] ) {
        	WarningOut().fmt("Input file with version: '{}' is no compatible with the program version: '{}' \n",
                    input_version, version);
        }

        // should flow123d wait for pressing "Enter", when simulation is completed
        sys_info.pause_after_run = i_rec.val<bool>("pause_after_run");
        // read record with problem configuration
        Input::AbstractRecord i_problem = i_rec.val<AbstractRecord>("problem");

        if (i_problem.type() == HC_ExplicitSequential::get_input_type() ) {

            problem_ = new HC_ExplicitSequential(i_problem);

            // run simulation
            problem_->run_simulation();
        } else {
            THROW( ExcUnknownProblem() );
        }

    }

    this->after_run();
}




void Application::after_run() {
	if (sys_info.pause_after_run) {
        printf("\nPress <ENTER> for closing the window\n");
        getchar();
    }
}


void _transform_profiler_data (const string &json_filepath, const string &output_file_suffix, const string &formatter) {
	namespace py = pybind11;

    if (json_filepath == "") return;

    // grab module and function by importing module profiler_formatter_module.py
    auto python_module = PythonLoader::load_module_by_name ("py123d.profiler.profiler_formatter_module");
    //
    // def convert (json_location, output_file, formatter):
    //
    auto convert_method = python_module.attr("convert");
    // execute method with arguments
    convert_method(json_filepath, (json_filepath + output_file_suffix), formatter);

}





Application::~Application() {
	if (problem_) delete problem_;

    if (use_profiler) {
    	// TODO: make a static output method that does nothing if the instance does not exist yet.
    	string profiler_json;
        if (petsc_initialized) {
            // log profiler data to this stream
            profiler_json = Profiler::instance()->output(PETSC_COMM_WORLD, profiler_path);
        } else {
        	profiler_json = Profiler::instance()->output(profiler_path);
        }

        // call python script which transforms json file at given location
        // Profiler::instance()->transform_profiler_data (".csv", "CSVFormatter");
        _transform_profiler_data (profiler_json, ".txt", "SimpleTableFormatter2");

        // finally uninitialize
        Profiler::uninitialize();
    }

    // TODO: have context manager classes for petsc and permon initialization
    // create a local variable in Application::run or so
	permon_finalize();
	petcs_finalize();
}




