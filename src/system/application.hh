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
 * @file    application_base.hh
 * @brief   
 */

#ifndef APPLICATION_BASE_HH_
#define APPLICATION_BASE_HH_

#include "system/application.hh"
//#include <mpi.hh>
//#include "system/system.hh"
//#include "system/sys_profiler.hh"
#include "system/python_loader.hh"
#include "coupling/hc_explicit_sequential.hh"
#include "coupling/balance.hh"
//#include "input/accessors.hh"
//#include "input/reader_to_storage.hh"
//#include "input/reader_internal_base.hh"
#include "system/armadillo_tools.hh"
//
#include <iostream>
#include <fstream>
#include <regex>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
//#include <boost/filesystem.hpp>


#include <string>
//#include <sstream>
#include <mpi.h>
//#include "global_defs.h"
//#include "system/file_path.hh"

#include <stdarg.h>                  // for va_list
#include <stdio.h>                   // for FILE

#include "config.h"                  // for FLOW123D_HAVE_PETSC
#include "petscsys.h"                // for PetscErrorCode
#include "system/exceptions.hh"      // for ExcStream, operator<<, EI, TYPED...

#include "input/type_output.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "input/reader_internal_base.hh"


//#ifdef FLOW123D_HAVE_PETSC
//#include "petsc.h"
//#endif

using namespace std;

// Exception that prints the signal number and name.
TYPEDEF_ERR_INFO( EI_Signal, int);
TYPEDEF_ERR_INFO( EI_SignalName, string);
DECLARE_EXCEPTION( ExcSignal, << "[ Signal " << EI_Signal::val << " (" << EI_SignalName::val << ") received! ]" );




/**
 * Base virtual class of Flow123D application.
 *
 * Contains base methods of application for initialization, run and finalization which is used in derived class.
 *
 * Usage:
 @code
 class Application : public ApplicationBase {
 public:
   Application(int argc, char ** argv); // constructor
 protected:
   virtual void run(); // implementation of pure virtual method
   virtual void after_run(); // overriding of parent method
 }
 @endcode
 *
 */
class Application {
public:
    TYPEDEF_ERR_INFO( EI_InputVersionStr, string);
    DECLARE_EXCEPTION( ExcVersionFormat,
            << "Wrong format of the version specification: "
            << EI_InputVersionStr::qval);
    DECLARE_INPUT_EXCEPTION( ExcUnknownProblem, << "Problem type not implemented.\n" );
    DECLARE_INPUT_EXCEPTION( ExcNoRunOption, << "No run option should be catched. Seeng this message is an error.\n" );

    /// Return codes of application
	static const int exit_success = 0;
    static const int exit_failure = 1;
    static const int exit_output = 0;	//return code if printout (text, JSON or LaTeX) is run

    static bool petsc_initialized;
    static bool permon_initialized;



    /// Root of the Input::Type tree. Description of whole input structure.
    static Input::Type::Record & get_input_type();


    /**
	 * Constructor
	 *
	 * Construction is done in init method. We need to call virtual methods during construction.
	 */
	Application();


    /**
     * Read main input file
     *
     * Returns accessor to the root Record.
     */
    Input::Record read_input();



    void init(int argc, char ** argv);
    /**
     * Run application.
     *
     * Read input and solve problem.
     */
    void run();


    /**
     * Displays program version and build info.
     * Pass version information to Profiler.
     *
     * TODO: Split these two functionalities.
     */
    void display_version();


    /// Destructor
    ~Application();


protected:

    /**
     * Check pause_after_run flag defined in input file.
     */
    void after_run();



    /**
     * Parse command line parameters.
     * @param[in] argc       command line argument count
     * @param[in] argv       command line arguments
     */
     void parse_cmd_line(const int argc, char ** argv);

protected:



	/**
	 * Read system parameters, open log.
	 */
	void system_init( MPI_Comm comm, const string &log_filename);


	/**
	 * Implement printf function for PETSc with support for redirecting.
	 */
#ifdef FLOW123D_HAVE_PETSC
	static PetscErrorCode petscvfprintf(FILE *fd, const char format[], va_list Argp);
#endif

	/**
	 * Initialize PETSC.
	 */
	void petsc_initialize(int argc, char ** argv);

	/**
	 * Finalize PETSC. If finalization failed return nonzero value.
	 */
	int petcs_finalize();

	/**
	 * Initialize PERMON.
	 */
	void permon_initialize(int argc, char ** argv);

	/**
	 * Finalize PERMON. If finalization failed return nonzero value.
	 */
	int permon_finalize();


    /**
     * Log file name argument - passed to system_init; "" means default, "\n" means no logging
     * TODO: move whole system_init into Application, use singleton for some runtime global options
     * for the Flow123d library.
     */
    string log_filename_;

    ///  Optional file name for output of PETSc parameters.
    ///  Has to be set in @p parse_cmd_line()
    string petsc_redirect_file_="";

    /// File handler for redirecting PETSc output
    static FILE *petsc_output_;

    /// Turn off signal handling useful to debug with valgrind.
    bool signal_handler_off_;


    /// Get version of program and other base data from rev_num.h and store them to map
    //Input::Type::RevNumData get_rev_num_data();

    /// Main Flow123d problem
    HC_ExplicitSequential *problem_;

    /// filename of main input file
    string main_input_filename_;

    //int passed_argc_;
    //char ** passed_argv_;

    /// Description of possible command line arguments.
    string program_arguments_desc_;

    /// If true, we do output of profiling information.
    bool use_profiler;

    /// location of the profiler report file
    string profiler_path;

    /// If true, preserves output of balance in YAML format.
    bool yaml_balance_output_;

    /// root input record
    Input::Record root_record;

};






#endif /* APPLICATION_BASE_HH_ */
