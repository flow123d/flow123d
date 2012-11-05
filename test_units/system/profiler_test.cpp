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
#include <ctime>

#include "system/sys_profiler.hh"

TEST(Profiler, basic_usage) {

    Profiler::initialize(MPI_COMM_WORLD);
    START_TIMER("test_tag");

    int j;
    for(int i=0;i<1000000;i++) j=2*i;

    END_TIMER("test_tag");
    Profiler::uninitialize();
}


TEST(Profiler, clock_timing) {
    Profiler::initialize(MPI_COMM_WORLD);
    unsigned int cycles = 1000000;

    clock_t t1,t2;

    t2=t1=clock();
    while (t1 == t2) { clock(); t2=clock();}
    cout << "Minimum timer resolution: " << t2-t1 << "per sec:" << CLOCKS_PER_SEC << endl;

    clock_t start;
    t1= clock();
    for(int i=0;i<cycles;i++) {
        start=clock(); start=clock();
    }
    t2=clock();
    cout << "ticks per 2*10^6 times clock(): " << t2-t1 << endl;

    start = clock();
    for(int i=0;i<cycles;i++) {
        START_TIMER("tag");
        END_TIMER("tag");
    }
    cout << "ticks per timerFrame: " << clock() - start << endl;
    Profiler::uninitialize();
}



TEST(Profiler, structure) {
    Profiler::initialize(MPI_COMM_WORLD);

    {
        START_TIMER("main");

        START_TIMER("sub1");
           START_TIMER("cross");
        END_TIMER("sub1");

        START_TIMER("sub2");
            END_TIMER("cross");
            START_TIMER("sub_sub");
                START_TIMER("sub1");
                END_TIMER("sub1");
            END_TIMER("sub_sub");
        END_TIMER("sub2");


    }
    Profiler::uninitialize();

}
