/*
 * profiler_test.cpp
 *
 *  Created on: Jun 16, 2012
 *      Author: jb
 */

#include <ctime>
#include <cstdlib>
#include <sstream>



#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>



#include "system/system.hh"
#include "system/sys_profiler.hh"

#ifdef DEBUG_PROFILER

/*************************
 * We test collisions of hash function on strings with max length 13.
 */

unsigned int random_string(char *str){
    unsigned int length = rand()%12+1;  // random string length from 1 up to 13 characters
    unsigned int i;
    for(i=0; i < length; i++) {
        str[i] = rand()%(128-32)+32; // random character from 32 (space) till the and of ASCII
    }
    str[i]=0;
    //printf("str: '%s'\n", str);
    return length;
}

TEST(Profiler, str_hash) {
    EXPECT_EQ(0, str_hash(""));
    EXPECT_EQ(65, str_hash("A"));
    EXPECT_EQ(6597, str_hash(" A"));

    srand ( time(NULL) );
    char a[16];
    char b[16];

// random test for hash collision
    unsigned int n_pairs=100;
    for(unsigned int i=0; i<n_pairs; i++) {
        random_string(a);
        random_string(b);
        if (string(a) != string(b) )
            EXPECT_NE( str_hash(a) , str_hash(b) );
    }
}



TEST(Profiler, CodePoint) {
    CodePoint cp = CODE_POINT("my_tag");    unsigned int line_save = __LINE__;
    EXPECT_EQ("my_tag", cp.tag_);
    EXPECT_EQ( str_hash("my_tag"), cp.hash_);
    EXPECT_EQ( line_save, cp.line_ );
    //EXPECT_EQ(string("(profiler_test.cpp, xxx(), 130)"), string(CODE_POINT) );
}




/**********************************************
 * Tests of Profiler.
 */

// wait given amount of time (in ms) and return it in ms
double wait( double time) {
//    cout << "wait function\n" <<endl;
    clock_t t1,t2;
    clock_t int_time = time /1000.0 * CLOCKS_PER_SEC;
    t2=t1=clock();
    
    while (t1 + int_time > t2) { t2=clock();}
    double time_in_ms = 1000.0 * (t2-t1) / CLOCKS_PER_SEC;
    cout << "time to wait: " << time << " actual: " << time_in_ms  << endl;
    
    return time_in_ms;
}

// return smallest amount of time resoluted by clock() function
double clock_resolution() {
//    cout << "wait function\n" <<endl;
    clock_t t1,t2;

    t2=t1=clock();
    while (t1 == t2) { clock(); t2=clock();}
    double min_time_step = 1000.0 * (t2-t1) / CLOCKS_PER_SEC;
    cout << "min wait time: " << min_time_step << endl;
    return min_time_step;
}


#define AT    string(Profiler::instance()->actual_tag())
#define ACT    Profiler::instance()->actual_cumulative_time()
#define AC    Profiler::instance()->actual_count()

TEST(Profiler, one_timer) {

    Profiler::initialize();
    double wait_time = clock_resolution();
    double total=0;
    { // uninitialize can not be in the same block as the START_TIMER
    START_TIMER("test_tag");
    EXPECT_EQ( 1, AC);
    total += wait(wait_time);
    END_TIMER("test_tag");

    START_TIMER("test_tag");
    EXPECT_EQ( total, ACT);
    EXPECT_EQ( 2, AC);
    total += wait(wait_time);
    total += wait(wait_time);
    END_TIMER("test_tag");

    START_TIMER("test_tag");
    EXPECT_EQ( total, ACT);
    EXPECT_EQ( 3, AC);

    }

    // test add_call
    {
        START_TIMER("add_call");
        ADD_CALLS(1000);
        EXPECT_EQ(1000, AC);
    }

    Profiler::uninitialize();

}




/*
  // This is efficiency test of START_TIMER macro.
  // It will pass only with optimalized build (not debug).


TEST(Profiler, efficiency) {
    Profiler::initialize();
    unsigned int cycles = 500000;

    double min_time_period = clock_resolution();
    cout << "Minimum timer resolution: " << min_time_period << " ms" << endl;
    EXPECT_LT( min_time_period, 20 ); // resolution better then 20ms

    double t_clock, t_timer;
    // measure the clock() function
    {
        START_TIMER("clock");
        for(int i=0;i<cycles;i++) {
            clock(); clock();  // two measurements are necessary for one frame
        }
    }

    // measure Timer overhead
    {
        START_TIMER("timer");
        for(int i=0;i<cycles;i++) {
            START_TIMER("tag");
            END_TIMER("tag");
        }

    }

    // get both cumulative times
    { START_TIMER("clock");
      t_clock = ACT;
    }

    { START_TIMER("timer");
      t_timer = ACT;
    }

    cout << "ticks per 2*10^6 times clock(): " << t_clock << endl;
    cout << "ticks per timer frame: " << t_timer << endl;
    EXPECT_LT( t_timer , 2*t_clock); // we allow timer to run two times slower compared to clock() function itself.

    Profiler::uninitialize();
}

*/


TEST(Profiler, structure) {
    Profiler::initialize();

    {
        START_TIMER("main");
        EXPECT_EQ("main", AT );

        START_TIMER("sub1");
           EXPECT_EQ("sub1", AT);
           START_TIMER("cross");
           EXPECT_EQ("cross", AT);

        END_TIMER("sub1");
        EXPECT_EQ("main", AT );


        START_TIMER("sub2");
            END_TIMER("cross");
            START_TIMER("sub_sub");
                START_TIMER("sub1");
                END_TIMER("sub1");
            END_TIMER("sub_sub");
        END_TIMER("sub2");

        START_TIMER("sub1");
        END_TIMER("sub1");
    }

    std::stringstream sout;
    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::instance()->output(MPI_COMM_WORLD, sout);

    EXPECT_TRUE( sout.str().find("Whole Program   0") );
    EXPECT_TRUE( sout.str().find("  sub1          2") );

    Profiler::uninitialize();

}

#else // DEBUG_PROFILER


TEST(Profiler, test_calls_only) {
    Profiler::initialize();
    START_TIMER("sub1");
    END_TIMER("sub1");
    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();

}


#endif // DEBUG_PROFILER


