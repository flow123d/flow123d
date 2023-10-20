/*
 * fem_tools_bench.cpp
 *
 *  Created on: Oct 18, 2023
 *      Author: David Flanderka
 *
 *  Speed tests of determinant function
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <armadillo>
#include "fem/fem_tools.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"


class FemToolsTest : public testing::Test {
public:
	FemToolsTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/fem";
        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::set_memory_monitoring(false, false);
    }

    ~FemToolsTest()
    {
        Profiler::uninitialize();
    }

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};



/// Check correct implementation of 'determinant()' function
//TEST(FemToolsDevelopTest, determinant) {
//    arma::mat::fixed<1,1> mat11 = {2};
//    arma::mat::fixed<2,2> mat22 = { {2, 3}, {4, 5} };
//    arma::mat::fixed<3,3> mat33 = { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} };
//
//    EXPECT_DOUBLE_EQ( det(mat11), determinant(mat11) );
//    EXPECT_DOUBLE_EQ( det(mat22), determinant(mat22) );
//    EXPECT_DOUBLE_EQ( det(mat33), determinant(mat33) );
//}


/// Benchmark test. Compare 'determinant()' and 'arma::det()' function. Matrix size: 3x3
TEST_F(FemToolsTest, speed_det_mat33) {
    static const uint N_RUNS = 2.5e7;

    std::vector< arma::mat::fixed<3,3> > mat_vec = {
            { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} },
            { {2, 4, 5}, {3, 2, 1}, {0, 6, 4} },
            { {1, 4, 2}, {5, 1, 3}, {2, 3, 4} },
            { {9, 7, 5}, {2, 4, 6}, {1, 3, 8} }
    };

    uint vec_size = mat_vec.size();
    std::vector< double > result_own(vec_size);
    std::vector< double > result_arma(vec_size);

    START_TIMER("own_implementation");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<mat_vec.size(); ++j) result_own[j] = determinant( mat_vec[j] );
    END_TIMER("own_implementation");

    START_TIMER("arma_implementation");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<mat_vec.size(); ++j) result_own[j] = det( mat_vec[j] );
    END_TIMER("arma_implementation");

    this->profiler_output("fem_tools");
}

