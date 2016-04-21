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

    LOG(_warning) << "Test of logger\n - type: warning";
}
