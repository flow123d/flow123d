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
#include "petscvec.h"  
#include "petscsys.h"  

#ifdef FLOW123D_DEBUG_PROFILER

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
    EXPECT_EQ(0, str_hash("", PROFILER_HASH_DEFAULT));
    EXPECT_EQ(65, str_hash("A", PROFILER_HASH_DEFAULT));
    EXPECT_EQ(6597, str_hash(" A", PROFILER_HASH_DEFAULT));

    srand ( time(NULL) );
    char a[16];
    char b[16];

// random test for hash collision
    unsigned int n_pairs=100;
    for(unsigned int i=0; i<n_pairs; i++) {
        random_string(a);
        random_string(b);
        if (string(a) != string(b) )
            EXPECT_NE( str_hash(a, PROFILER_HASH_DEFAULT) , str_hash(b, PROFILER_HASH_DEFAULT) );
    }
}



TEST(Profiler, CodePoint) {
    CodePoint cp = CODE_POINT("my_tag");    unsigned int line_save = __LINE__;
    EXPECT_EQ("my_tag", cp.tag_);
    EXPECT_EQ( str_hash("my_tag", PROFILER_HASH_DEFAULT), cp.hash_);
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

// wait given amount of time (in sec) and return it in sec
double wait_sec( double time) {
	TimePoint t1, t2;

	t2 = t1 = TimePoint();
	while ((t2-t1) < time)  {
		t2 = TimePoint();
//		cout << "difference: " << (t2-t1) << endl;
	}

	return (t2-t1);
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
#define ACT   Profiler::instance()->actual_cumulative_time()
#define AC    Profiler::instance()->actual_count()


TEST(Profiler, one_timer) {

//	TimePoint t0 = TimePoint();
//	wait(1000);
//	TimePoint t1 = TimePoint();
//
//	cout << "CURRENT " << t1-t0 << endl;

	const double TIMER_RESOLUTION = Profiler::get_resolution();
	const double DELTA = TIMER_RESOLUTION*100;
	double total=0;
    Profiler::initialize();


    { // uninitialize can not be in the same block as the START_TIMER


    START_TIMER("test_tag");
    	// test that number of calls of current timer is
	    EXPECT_EQ( 1, AC);

	    // wait a TIMER_RESOLUTION time
		total += wait_sec(TIMER_RESOLUTION);
    END_TIMER("test_tag");



    START_TIMER("test_tag");
    	// test that number of calls of current timer is
		EXPECT_EQ( 2, AC);

		// test whether difference between measured time and total time is within TIMER_RESOLUTION
		EXPECT_LE( abs(ACT-total), DELTA);
		cout << "difference: " << abs(total-ACT) << ", tolerance: " << DELTA << endl;

		// wait a TIMER_RESOLUTION time
		total += wait_sec (TIMER_RESOLUTION);
		total += wait_sec (TIMER_RESOLUTION);

    END_TIMER("test_tag");

//    for (int i = 0; i < 100; i++) {
//        START_TIMER("test_tag");
//            EXPECT_LE( abs(ACT-total), DELTA);
//            cout << i+1 <<". difference: " << abs(total-ACT) << ", tolerance: " << DELTA << endl;
//            total += wait_sec (TIMER_RESOLUTION);
//    	END_TIMER("test_tag");
//    }

    START_TIMER("test_tag");
    	EXPECT_EQ( 3, AC);
		EXPECT_LE( abs(ACT-total), DELTA);
		cout << "difference: " << abs(total-ACT) << ", tolerance: " << DELTA << endl;
    }



    // test add_call
    {
        START_TIMER("add_call");
        ADD_CALLS(1000);
        EXPECT_EQ(1000, AC);
    }

    // test absolute time
    {
        START_TIMER("one_second");
        wait_sec(1);
    }
    std::stringstream sout;
    Profiler::instance()->output(MPI_COMM_WORLD, sout);
    Profiler::instance()->output(MPI_COMM_WORLD, cout);

    //EXPECT_NE( sout.str().find("\"tag\": \"Whole Program\""), string::npos );

    Profiler::uninitialize();
}



TEST(Profiler, absolute_time) {
    Profiler::initialize();

    // test absolute time
    {
        START_TIMER("one_second");
        wait_sec(1);
    }
    std::stringstream sout;
    Profiler::instance()->output(MPI_COMM_WORLD, sout);
    Profiler::instance()->output(MPI_COMM_WORLD, cout);

    EXPECT_NE( sout.str().find("cumul-time-min\": \"1.00"), string::npos );
    EXPECT_NE( sout.str().find("cumul-time-max\": \"1.00"), string::npos );

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

    // when using find, we need to compare result to string::npos value (which indicates not found)
    EXPECT_NE( sout.str().find("\"tag\": \"Whole Program\""), string::npos );
    EXPECT_NE( sout.str().find("\"tag\": \"sub1\""), string::npos );

    Profiler::uninitialize();

}



template <class T>
int alloc_and_dealloc(int size){
    T* t = new T[size];
    delete [] t;
    return size * sizeof(T);
}

TEST(Profiler, memory_add_calls) {
    int allocated = 0;
    Profiler::initialize();
    allocated += alloc_and_dealloc<int>(250);
    {
        for (int i = 0; i < 10000; i++) {
            START_TIMER("A");
            allocated += alloc_and_dealloc<int>(25);
            END_TIMER("A");
        }
    }
    EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
    Profiler::uninitialize();
}

TEST(Profiler, memory_profiler) {
    const int ARR_SIZE = 100;
    const int LOOP_CNT = 100;
    Profiler::initialize();

    {
        START_TIMER("memory-profiler-int");
        // alloc and dealloc array of int
        for (int i = 0; i < LOOP_CNT; i++) alloc_and_dealloc<int>(ARR_SIZE);
        // test that we deallocated all allocated space
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), Profiler::instance()->actual_memory_dealloc());
        // test that allocated space is correct size
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), ARR_SIZE * LOOP_CNT * sizeof(int));
        END_TIMER("memory-profiler-int");


        START_TIMER("memory-profiler-double");
        // alloc and dealloc array of float
        for (int i = 0; i < LOOP_CNT; i++) alloc_and_dealloc<double>(ARR_SIZE);
        // test that we deallocated all allocated space
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), Profiler::instance()->actual_memory_dealloc());
        // test that allocated space is correct size
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), ARR_SIZE * LOOP_CNT * sizeof(double) + ARR_SIZE * LOOP_CNT * sizeof(int));
        END_TIMER("memory-profiler-double");


        START_TIMER("memory-profiler-simple");
        // alloc and dealloc array of float
        for (int i = 0; i < LOOP_CNT; i++) {
            int * j = new int;
            delete j;
        }
        // test that we deallocated all allocated space
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), Profiler::instance()->actual_memory_dealloc());
        // test that allocated space is correct size
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), LOOP_CNT * sizeof(int) + ARR_SIZE * LOOP_CNT * sizeof(double) + ARR_SIZE * LOOP_CNT * sizeof(int));
        END_TIMER("memory-profiler-simple");
    }

    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

TEST(Profiler, memory_propagation){
    const int SIZE = 25;
    int allocated = 0;
    
    Profiler::initialize();
    {
        allocated += alloc_and_dealloc<int>(SIZE);
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
        
        START_TIMER("A");
            allocated += alloc_and_dealloc<int>(10 * SIZE);
            EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
            
            START_TIMER("B");
                allocated += alloc_and_dealloc<int>(100 * SIZE);
                EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
                
                START_TIMER("C");
                    EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
                END_TIMER("C");
                
                EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
            END_TIMER("B");
            
            allocated += alloc_and_dealloc<int>(10 * SIZE);
            EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
            
            START_TIMER("D");
                EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
            END_TIMER("D");
            
            EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
        END_TIMER("A");
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
    }
    EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
    EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), Profiler::instance()->actual_memory_dealloc());
    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

TEST(Profiler, petsc_memory_monitor) {
    PetscInt size = 10000;
    Profiler::initialize();
        
    {
        START_TIMER("A");
            Vec tmp_vector;
            VecCreateSeq(MPI_COMM_WORLD, size, &tmp_vector);
            VecDestroy(&tmp_vector);
            
            START_TIMER("C");
            END_TIMER("C");
        END_TIMER("A");
        
        START_TIMER("B");
            Vec tmp_vector1, tmp_vector2;
            VecCreateSeq(MPI_COMM_WORLD, size, &tmp_vector1);
            VecCreateSeq(MPI_COMM_WORLD, size, &tmp_vector2);
            VecDestroy(&tmp_vector1);
            VecDestroy(&tmp_vector2);
        END_TIMER("B");
    }
    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

TEST(Profiler, multiple_instances) {
    int allocated = 0;
    for (int i = 0; i < 5; i++) {
        allocated = 0;
        Profiler::initialize();
        {
            allocated += alloc_and_dealloc<int>(25);
        }
        EXPECT_EQ(Profiler::instance()->actual_memory_alloc(), allocated);
        Profiler::uninitialize();
    }
}


#else // FLOW123D_DEBUG_PROFILER


TEST(Profiler, test_calls_only) {
    Profiler::initialize();
    START_TIMER("sub1");
    END_TIMER("sub1");
    Profiler::instance()->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();

}


#endif // FLOW123D_DEBUG_PROFILER
