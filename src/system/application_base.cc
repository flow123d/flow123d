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

#include "system/application_base.hh"
#include "system/sys_profiler.hh"
#include "system/logger_options.hh"
#include "system/armadillo_tools.hh"
#include "system/file_path.hh"
#include "system/system.hh"
#include <signal.h>

#ifdef FLOW123D_HAVE_PETSC
//#include <petsc.h>
#include <petscsys.h>
#include <petsc/private/petscimpl.h> /* to gain access to the private PetscVFPrintf */
#endif

#include <string.h>                                    // for strsignal
#include <boost/exception/detail/error_info_impl.hpp>  // for error_info
#include <iostream>                                    // for cout
#include <sstream>                                     // for operator<<, endl
#include "mpi.h"                                       // for MPI_Comm_size
#include "petscerror.h"                                // for CHKERRQ, Petsc...
#include "system/exc_common.hh"                        // for ExcAssertMsg
#include "system/global_defs.h"                        // for OLD_ASSERT, msg
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


ApplicationBase::ApplicationBase()
: log_filename_(""),
  signal_handler_off_(false)
{ }

bool ApplicationBase::petsc_initialized = false;


void ApplicationBase::system_init( MPI_Comm comm, const string &log_filename ) {
    int ierr;

    sys_info.comm=comm;


    //Xio::init(); //Initialize XIO library

    // TODO : otevrit docasne log file jeste pred ctenim vstupu (kvuli zachyceni chyb), po nacteni dokoncit
    // inicializaci systemu

    ierr=MPI_Comm_rank(comm, &(sys_info.my_proc));
    ierr+=MPI_Comm_size(comm, &(sys_info.n_proc));
    LoggerOptions::get_instance().setup_mpi(comm);
    OLD_ASSERT( ierr == MPI_SUCCESS,"MPI not initialized.\n");

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


FILE *ApplicationBase::petsc_output_ =NULL;

#ifdef FLOW123D_HAVE_PETSC
PetscErrorCode ApplicationBase::petscvfprintf(FILE *fd, const char format[], va_list Argp) {
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


void ApplicationBase::petsc_initialize(int argc, char ** argv) {
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



int ApplicationBase::petcs_finalize() {
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


void ApplicationBase::init(int argc, char ** argv) {
    // parse our own command line arguments, leave others for PETSc
	this->parse_cmd_line(argc, argv);
    Profiler::instance();

    armadillo_setup(); // set catching armadillo exceptions and reporting stacktrace

	this->petsc_initialize(argc, argv);
	petsc_initialized = true;

    this->system_init(PETSC_COMM_WORLD, log_filename_); // Petsc, open log, read ini file
}


ApplicationBase::~ApplicationBase() {
	//if (sys_info.log) xfclose(sys_info.log);
	petcs_finalize();
}
