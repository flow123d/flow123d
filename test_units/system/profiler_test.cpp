/*
 * profiler_test.cpp
 *
 *  Created on: Jun 16, 2012
 *      Author: jb
 */


#define DEBUG
#define DEBUG_PROFILER

#define TEST_USE_MPI
#include <gtest_mpi.hh>

#include "system/sys_profiler.hh"

TEST(Profiler, basic_usage) {

    Profiler::initialize(MPI_COMM_WORLD);
    START_TIMER("test_tag");

    int j;
    for(int i=0;i<1000000;i++) j=2*i;

    END_TIMER("test_tag");
    Profiler::uninitialize();
}
