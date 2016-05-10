/*
 * logger_test.cpp
 *
 *  Created on: May 10, 2012
 *      Author: jb
 */


#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/logger.hh"
#include "system/sys_profiler.hh"


TEST(Logger, CompleteTest) {
	Profiler::initialize();
	LoggerOptions::get_instance().setup_mpi(MPI_COMM_WORLD);
	LoggerOptions::get_instance().set_log_file("flow123d");

	MessageOut() << "Test of logger \n... next line" << "\n" << "... next line" << std::endl;

	MessageOut().every_proc() << "Output with flag every process." << std::endl;

	unsigned int mesh_size = 150;
	double start_time = 0.5;
	WarningOut() << "Start of simulation at time " << start_time << ", mesh has " << mesh_size
    		<< " elements." << std::endl << "... next line in separate message" << std::endl;
}
