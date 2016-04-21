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
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

/**
 * Fixture class for testing Profiler protected members
 */
class ProfilerTest;

#define __UNIT_TEST__
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "petscvec.h"
#include "petscsys.h"



class ProfilerTest: public testing::Test {
    public:
        void test_str_hash();
        void test_code_point();
        void test_one_timer();
        void test_absolute_time();
        void test_structure();
        void test_memory_profiler();
        void test_petsc_memory();
        void test_memory_propagation();
        void test_petsc_memory_monitor();
        void test_multiple_instances();
        void test_propagate_values();
        // void test_inconsistent_tree();
};


#ifdef FLOW123D_DEBUG_PROFILER


#define PI       Profiler::instance()
#define AN       PI->timers_[PI->actual_node]
#define ATN      string(AN.tag())
#define ACT      AN.cumulative_time()
#define ACC      AN.call_count
#define MALLOC   AN.total_allocated_
#define DEALOC   AN.total_deallocated_


// ---------------- //
// Helper functions //
// ---------------- //


// populate given str with random string and return its length <1,12>
unsigned int random_string(char *str){
    unsigned int length = rand()%12+1;  // random string length from 1 up to 13 characters
    unsigned int i;
    for(i=0; i < length; i++) {
        str[i] = rand()%(128-32)+32; // random character from 32 (space) till the and of ASCII
    }
    str[i]=0;
    return length;
}

// simple function for allocating and deallocating array <T> of given length
template <class T>
int alloc_and_dealloc(int size){
    T* t = new T[size];
    delete [] t;
    return size * sizeof(T);
}

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
//      cout << "difference: " << (t2-t1) << endl;
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


// ----------------------------- //
// Profiler test implementations //
// ----------------------------- //


// We test collisions of hash function on strings with max length 13.
TEST_F(ProfilerTest, test_str_hash) {test_str_hash();}
void ProfilerTest::test_str_hash() {
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

// testing CodePoint contains correct value
TEST_F(ProfilerTest, test_code_point) {test_code_point();}
void ProfilerTest::test_code_point() {
    CodePoint cp = CODE_POINT("my_tag");    unsigned int line_save = __LINE__;
    EXPECT_EQ("my_tag", cp.tag_);
    EXPECT_EQ( str_hash("my_tag", PROFILER_HASH_DEFAULT), cp.hash_);
    EXPECT_EQ( line_save, cp.line_ );
    //EXPECT_EQ(string("(profiler_test.cpp, xxx(), 130)"), string(CODE_POINT) );
}

// testing profiler precision up to 2 decimal places relative to TIMER_RESOLUTION
TEST_F(ProfilerTest, test_one_timer) {test_one_timer();}
 void ProfilerTest::test_one_timer() {
    const double TIMER_RESOLUTION = Profiler::get_resolution();
    const double DELTA = TIMER_RESOLUTION*1000;
    double total=0;
    Profiler::initialize();

    { // uninitialize can not be in the same block as the START_TIMER


    START_TIMER("test_tag");
        // test that number of calls of current timer is
        EXPECT_EQ( 1, ACC);

        // wait a TIMER_RESOLUTION time
        total += wait_sec(TIMER_RESOLUTION);
    END_TIMER("test_tag");



    START_TIMER("test_tag");
        // test that number of calls of current timer is
        EXPECT_EQ( 2, ACC);

        // test whether difference between measured time and total time is within TIMER_RESOLUTION
        EXPECT_LE( abs(ACT-total), DELTA);
        cout << "difference: " << abs(total-ACT) << ", tolerance: " << DELTA << endl;

        // wait a TIMER_RESOLUTION time
        total += wait_sec (TIMER_RESOLUTION);
        total += wait_sec (TIMER_RESOLUTION);

    END_TIMER("test_tag");

    START_TIMER("test_tag");
        EXPECT_EQ( 3, ACC);
        EXPECT_LE( abs(ACT-total), DELTA);
        cout << "difference: " << abs(total-ACT) << ", tolerance: " << DELTA << endl;
    }

    // test add_call
    {
        START_TIMER("add_call");
        ADD_CALLS(1000);
        EXPECT_EQ(1000, ACC);
    }

    // test absolute time
    {
        START_TIMER("one_second");
        wait_sec(1);
    }
    std::stringstream sout;
    PI->output(MPI_COMM_WORLD, sout);
    PI->output(MPI_COMM_WORLD, cout);

    //EXPECT_NE( sout.str().find("\"tag\": \"Whole Program\""), string::npos );

    Profiler::uninitialize();
}

// testing precision when waiting 1 sec up to 2 decimal places
TEST_F(ProfilerTest, test_absolute_time) {test_absolute_time();}
void ProfilerTest::test_absolute_time() {
    Profiler::initialize();

    // test absolute time
    {
        START_TIMER("one_second");
        wait_sec(1);
    }
    std::stringstream sout;
    PI->output(MPI_COMM_WORLD, sout);
    PI->output(MPI_COMM_WORLD, cout);
    
    // try to transform profiler data using python
    PI->output(MPI_COMM_WORLD);
    PI->transform_profiler_data (".txt", "SimpleTableFormatter");
    
    int ierr, mpi_rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    EXPECT_EQ( ierr, 0 );
    
    // 0 processor will have valid profiler report
    // other processors should have empty string only
    if (mpi_rank == 0) {
        // test timer resolution, requiring atleast 2 digit places
        EXPECT_NE( sout.str().find("cumul-time-min\": \"1.00"), string::npos );
        EXPECT_NE( sout.str().find("cumul-time-max\": \"1.00"), string::npos );
    } else {
        EXPECT_TRUE( sout.str().empty() );
    }
    
    

    Profiler::uninitialize();
}

// testing correct report generation
TEST_F(ProfilerTest, test_structure) {test_structure();}
void ProfilerTest::test_structure() {
    Profiler::initialize();

    {
        START_TIMER("main");
        EXPECT_EQ("main", ATN);

        START_TIMER("sub1");
           EXPECT_EQ("sub1", ATN);
           START_TIMER("cross");
           EXPECT_EQ("cross", ATN);

        END_TIMER("sub1");
        EXPECT_EQ("main", ATN );


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
    PI->output(MPI_COMM_WORLD, cout);
    PI->output(MPI_COMM_WORLD, sout);


    int ierr, mpi_rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    EXPECT_EQ( ierr, 0 );
    
    // 0 processor will have valid profiler report
    // other processors should have empty string only
    if (mpi_rank == 0) {
        // when using find, we need to compare result to string::npos value (which indicates not found)
        EXPECT_NE( sout.str().find("\"tag\": \"Whole Program\""), string::npos );
        EXPECT_NE( sout.str().find("\"tag\": \"sub1\""), string::npos );
    } else {
        EXPECT_TRUE( sout.str().empty() );
    }
    
    Profiler::uninitialize();

}

// testing memory alloc and dealloc in separated time-frames
TEST_F(ProfilerTest, test_memory_profiler) {test_memory_profiler();}
void ProfilerTest::test_memory_profiler() {
    const int ARR_SIZE = 1000;
    const int LOOP_CNT = 1000;
    Profiler::initialize();

    {
        START_TIMER("memory-profiler-int");
        // alloc and dealloc array of int
        for (int i = 0; i < LOOP_CNT; i++) alloc_and_dealloc<int>(ARR_SIZE);
        // test that we deallocated all allocated space
        EXPECT_EQ(MALLOC, DEALOC);
        // test that allocated space is correct size
        EXPECT_EQ(MALLOC, ARR_SIZE * LOOP_CNT * sizeof(int));
        END_TIMER("memory-profiler-int");


        START_TIMER("memory-profiler-double");
        // alloc and dealloc array of float
        for (int i = 0; i < LOOP_CNT; i++) alloc_and_dealloc<double>(ARR_SIZE);
        // test that we deallocated all allocated space
        EXPECT_EQ(MALLOC, DEALOC);
        // test that allocated space is correct size
        EXPECT_EQ(MALLOC, ARR_SIZE * LOOP_CNT * sizeof(double));
        END_TIMER("memory-profiler-double");


        START_TIMER("memory-profiler-simple");
        // alloc and dealloc array of float
        for (int i = 0; i < LOOP_CNT; i++) {
            int * j = new int;
            delete j;
        }
        // test that we deallocated all allocated space
        EXPECT_EQ(MALLOC, DEALOC);
        // test that allocated space is correct size
        EXPECT_EQ(MALLOC, LOOP_CNT * sizeof(int));
        END_TIMER("memory-profiler-simple");
    }

    PI->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

//testing simple petsc memory difference when manipulating with large data
TEST_F(ProfilerTest, test_petsc_memory) {test_petsc_memory();}
void ProfilerTest::test_petsc_memory() {
    int ierr, mpi_rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    EXPECT_EQ( ierr, 0 );
    
    Profiler::initialize(); {
        PetscLogDouble mem;
        START_TIMER("A");
            PetscInt size = 100*1000;
            PetscScalar value = 0.1;
            Vec tmp_vector;
            VecCreateSeq(PETSC_COMM_SELF, size, &tmp_vector);
            VecSet(tmp_vector, value);
            // VecSetRandom(tmp_vector, NULL);
        END_TIMER("A");
        
        START_TIMER("A");
            // allocated memory MUST be greater or equal to size * size of double
            EXPECT_GE(AN.petsc_memory_difference, size*sizeof(double));
        END_TIMER("A");
        
        START_TIMER("B");
            PetscScalar sum;
            VecSum(tmp_vector, &sum);
        END_TIMER("B");
        
        START_TIMER("C");
            VecDestroy(&tmp_vector);
        END_TIMER("C");
        
        START_TIMER("C");
            // since we are destroying vector, we expect to see negative memory difference
            EXPECT_LE(AN.petsc_memory_difference, 0);
        END_TIMER("C");
    }
    PI->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

//testing memory alloc and dealloc propagation in nested time-frames
TEST_F(ProfilerTest, test_memory_propagation) {test_memory_propagation();}
void ProfilerTest::test_memory_propagation(){
    const int SIZE = 25;
    int allocated_whole = 0;
    int allocated_A = 0;
    int allocated_B = 0;
    int allocated_C = 0;
    int allocated_D = 0;
    
    Profiler::initialize();
    {
        allocated_whole = MALLOC;
        allocated_whole += alloc_and_dealloc<int>(SIZE);
        EXPECT_EQ(MALLOC, allocated_whole);
        
        START_TIMER("A");
            allocated_A += alloc_and_dealloc<int>(10 * SIZE);
            EXPECT_EQ(MALLOC, allocated_A);
            
            START_TIMER("B");
                allocated_B += alloc_and_dealloc<int>(100 * SIZE);
                
                START_TIMER("C");
                    EXPECT_EQ(MALLOC, allocated_C);
                END_TIMER("C");
                allocated_B += allocated_C;
                
            END_TIMER("B");
            allocated_A += allocated_B;
            
            allocated_A += alloc_and_dealloc<int>(10 * SIZE);
            
            for(int i = 0; i < 5; i++) {
                START_TIMER("D");
                    allocated_D += alloc_and_dealloc<int>(1 * SIZE);
                END_TIMER("D");
                START_TIMER("D");
                    allocated_D += alloc_and_dealloc<int>(1 * SIZE);
                END_TIMER("D");
            }
            allocated_A += allocated_D;
            
            
        END_TIMER("A");
        allocated_whole += allocated_A;
    }
    PI->propagate_timers();
    EXPECT_EQ(MALLOC, allocated_whole);
    EXPECT_EQ(MALLOC, DEALOC);
    Profiler::uninitialize();
}

// testing petsc memory working properly
TEST_F(ProfilerTest, test_petsc_memory_monitor) {test_petsc_memory_monitor();}
void ProfilerTest::test_petsc_memory_monitor() {
    int ierr, mpi_rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    EXPECT_EQ( ierr, 0 );

    Profiler::initialize(); {
        PetscInt size = 10000;
        START_TIMER("A");
            Vec tmp_vector;
            VecCreateSeq(PETSC_COMM_SELF, size, &tmp_vector);
            VecDestroy(&tmp_vector);
            
            START_TIMER("C");
            END_TIMER("C");
        END_TIMER("A");
        
        START_TIMER("B");
            Vec tmp_vector1, tmp_vector2;
            VecCreateSeq(PETSC_COMM_SELF, size, &tmp_vector1);
            VecCreateSeq(PETSC_COMM_SELF, size, &tmp_vector2);
            VecDestroy(&tmp_vector1);
            VecDestroy(&tmp_vector2);
        END_TIMER("B");
    }
    PI->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

// testing multiple initialization and uninitialization of Profiler
TEST_F(ProfilerTest, test_multiple_instances) {test_multiple_instances();}
void ProfilerTest::test_multiple_instances() {
    int allocated = 0;
    for (int i = 0; i < 5; i++) {
        allocated = 0;
        Profiler::initialize();
        {
            allocated += alloc_and_dealloc<int>(25);
        }
        EXPECT_EQ(MALLOC, allocated);
        Profiler::uninitialize();
    }
}

// testing memory propagation with manual propagate_values call 
TEST_F(ProfilerTest, test_propagate_values) {test_propagate_values();}
void ProfilerTest::test_propagate_values() {
    int allocated = 0;
    Profiler::initialize(); {
            START_TIMER("A");
                START_TIMER("B");
                    START_TIMER("C");
                        allocated += alloc_and_dealloc<int>(25);
                    END_TIMER("C");
                END_TIMER("B");
                
                START_TIMER("D");
                END_TIMER("D");
                
                PI->propagate_timers();
                EXPECT_EQ(MALLOC, allocated);
            END_TIMER("A");
    }
    PI->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();
}

// optional test only for testing merging of inconsistent profiler trees
// TEST_F(ProfilerTest, test_inconsistent_tree) {test_inconsistent_tree();}
// void ProfilerTest::test_inconsistent_tree() {
//     // in order to fully test this case MPI consists of 2 processors at minimum
//     int mpi_rank;
//     std::stringstream sout;
//     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//     
//     Profiler::initialize();
//     if(mpi_rank == 0) {
//         START_TIMER("A");
//             START_TIMER("AA");
//             END_TIMER("AA");
//             
//             START_TIMER("BB");
//             END_TIMER("BB");
//         END_TIMER("A");
//     } else {
//         START_TIMER("A");
//             START_TIMER("AA");
//             END_TIMER("AA");
//         END_TIMER("A");
//         
//         START_TIMER("B");
//         END_TIMER("B");
//         
//         START_TIMER("C");
//         END_TIMER("C");
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//     PI->output(MPI_COMM_WORLD, sout);
//     
//     if (mpi_rank == 0) {
//         // tags BB B and C should not be part of report since they do not appear
//         // in all processors
//         EXPECT_EQ( sout.str().find("BB"), string::npos );
//         EXPECT_EQ( sout.str().find("B"),  string::npos );
//         EXPECT_EQ( sout.str().find("C"),  string::npos );
//         // tags A and AA are in all processors and must be in report
//         EXPECT_NE( sout.str().find("A"),  string::npos );
//         EXPECT_NE( sout.str().find("AA"), string::npos );
//     }
//     
//     Profiler::uninitialize();
// }



#else // FLOW123D_DEBUG_PROFILER

// testing non-fatal functioning of Profiler when debug is off
TEST(Profiler, test_calls_only) {
    Profiler::initialize();
    START_TIMER("sub1");
    END_TIMER("sub1");
    PI->output(MPI_COMM_WORLD, cout);
    Profiler::uninitialize();

}


#endif // FLOW123D_DEBUG_PROFILER
