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


TEST(Logger, String) {
	Profiler::initialize();

    LOG(_warning) << "Test of logger\n - type: warning" << "\n" << "Description" << std::endl;
}

TEST(Logger, CompleteTest) {
	Profiler::initialize();

	unsigned int mesh_size = 150;
	double start_time = 0.5;
    LOG(_message) << "Start of simulation at time " << start_time << ", mesh has " << mesh_size << " elements." << std::endl;
}
