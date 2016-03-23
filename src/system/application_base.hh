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


#include <string>
#include <sstream>
#include <mpi.h>

#include "global_defs.h"
#include "system/xio.h"
#include "system/file_path.hh"

#ifdef FLOW123D_HAVE_PETSC
#include "petsc.h"
#endif

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
class ApplicationBase {
public:

	/**
	 * Contains basic structure of application (initialization, run and finalization).
	 * Method is call after constructor and allows to call virtual methods.
	 */
	void init(int argc, char ** argv);

    /// Return codes of application
	static const int exit_success = 0;
    static const int exit_failure = 1;
    static const int exit_output = 0;	//return code if printout (text, JSON or LaTeX) is run

    static bool petsc_initialized;

protected:

	/**
	 * Constructor
	 *
	 * Construction is done in init method. We need to call virtual methods during construction.
	 */
	ApplicationBase(int argc,  char ** argv);

	/// Destructor
	virtual ~ApplicationBase();

	/**
	 * Run application.
	 *
	 * Method must be implemented in derived class.
	 */
	virtual void run() = 0;

	/**
	 * Read system parameters, open log.
	 */
	void system_init( MPI_Comm comm, const string &log_filename);

	/**
	 * Parse command line parameters before Flow123D initialization.
	 *
	 * Method can be override in derived class.
	 */
	virtual void parse_cmd_line(const int argc, char ** argv) {}

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
	 * Execute part of program after run of simulation.
	 *
	 * Method can be override in derived class.
	 */
	virtual void after_run() {}

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
};

#endif /* APPLICATION_BASE_HH_ */
