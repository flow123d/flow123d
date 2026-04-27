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


#include "coupling/application.hh"
#include <mpi.h>
//#include "system/system.hh"
//#include "system/sys_profiler.hh"
//#include "system/python_loader.hh"
//#include "coupling/hc_explicit_sequential.hh"
//#include "coupling/balance.hh"
//#include "input/accessors.hh"
//#include "input/reader_to_storage.hh"
//#include "input/reader_internal_base.hh"
//#include "system/armadillo_tools.hh"
//
#include <iostream>
//#include <fstream>
//#include <regex>
//#include <boost/program_options/parsers.hpp>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/options_description.hpp>
//#include <boost/filesystem.hpp>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds





void mpi_terminate() {
	int mpi_on;
	MPI_Initialized(&mpi_on);
	if (! mpi_on) return;

    // Test if all processes are in the exception.
    MPI_Request request;
    MPI_Ibarrier(MPI_COMM_WORLD, &request);
    std::this_thread::sleep_for(std::chrono::microseconds(10));
    int done;
    MPI_Status status;
    MPI_Test(&request, &done, &status);
    if (! done) {
        // Kill all if we can not synchronize.
        MPI_Abort( MPI_COMM_WORLD, Application::exit_failure);
    }
    // Peacefull end.
}


/**
 * Wrap application variable into separate function in order to
 * force correct call of destructors before main catch block which is necessary for
 * aborting all MPI processes.
 */
void application_run(int argc, char **argv) {
    Application app;
    app.init(argc, argv);
    app.run();
}

//=============================================================================

/**
 *  FUNCTION "MAIN"
 */
int main(int argc, char **argv) {
    try {
    	application_run(argc, argv);
    } catch (Application::ExcNoRunOption) {
        mpi_terminate();
        return Application::exit_success;
    } catch (std::exception & e) {
        _LOG( Logger::MsgType::error ).every_proc() << e.what();
        mpi_terminate();
        return Application::exit_failure;
    } catch (...) {
        _LOG( Logger::MsgType::error ).every_proc() << "Unknown exception" << endl;
        mpi_terminate();
        return Application::exit_failure;
    }

    // Say Goodbye
    return Application::exit_success;
}
