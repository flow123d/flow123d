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
#include "arma_expect.hh"

#include "fem/fem_tools.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"
#include "system/armor.hh"


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



/// Check correct implementation of 'determinant()' and 'inverse()' function
TEST(FemToolsDevelopTest, functionst) {
    arma::mat::fixed<1,1> mat11 = {2};
    arma::mat::fixed<2,2> mat22 = { {2, 3}, {4, 5} };
    arma::mat::fixed<3,3> mat33 = { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} };
    arma::mat::fixed<2,3> mat23 = { {1, 2, 3}, {4, 5, 6} };
    arma::mat::fixed<3,2> mat32 = { {1, 4}, {2, 5}, {3, 6} };

    /* Test of determinant function */
    EXPECT_DOUBLE_EQ( det(mat11), determinant(mat11) );
    EXPECT_DOUBLE_EQ( det(mat22), determinant(mat22) );
    EXPECT_DOUBLE_EQ( det(mat33), determinant(mat33) );

    /* Test of inverse function */
    // computing of inverse matrices
    arma::mat::fixed<1,1> inv11 = inverse(mat11);
    arma::mat::fixed<2,2> inv22 = inverse(mat22);
    arma::mat::fixed<3,3> inv33 = inverse(mat33);
    arma::mat::fixed<3,2> inv23 = inverse(mat23);
    arma::mat::fixed<2,3> inv32 = inverse(mat32);
    // expected values
    arma::mat::fixed<2,2> expect_22 = arma::eye(2,2);
    arma::mat::fixed<3,3> expect_33 = arma::eye(3,3);
    // matrix 1x1
    EXPECT_DOUBLE_EQ( mat11(0,0), 1 / inv11(0,0) );
    // matrix 2x2
    arma::mat::fixed<2,2> multi_22 = mat22 * inv22;
    EXPECT_ARMA_EQ(expect_22, multi_22);
    // matrix 3x3
    arma::mat::fixed<3,3> multi_33 = mat33 * inv33;
    EXPECT_ARMA_EQ(expect_33, multi_33);
    // matrix 2x3
    arma::mat::fixed<2,2> multi_23 = mat23 * inv23;
    EXPECT_ARMA_EQ(expect_22, multi_23);
    // matrix 3x2
    arma::mat::fixed<2,2> multi_32 = inv32 * mat32;
    EXPECT_ARMA_EQ(expect_22, multi_32);
}


/// Benchmark test. Compare 'determinant()' and 'arma::det()' function. Matrix size: 3x3
TEST_F(FemToolsTest, speed_test) {
    static const uint N_RUNS = 2.5e7;

    std::vector< arma::mat::fixed<3,3> > mat33_vec = {
            { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} },
            { {2, 4, 5}, {3, 2, 1}, {0, 6, 4} },
            { {1, 4, 2}, {5, 1, 3}, {2, 3, 4} },
            { {9, 7, 5}, {2, 4, 6}, {1, 3, 8} }
    };
    std::vector< arma::mat::fixed<2,3> > mat23_vec = {
            { {1, 2, 3}, {2, 4, 5} },
            { {2, 4, 5}, {3, 2, 1} },
            { {1, 4, 2}, {5, 1, 3} },
            { {9, 7, 5}, {2, 4, 6} }
    };

    uint vec_size = mat33_vec.size();
    std::vector< double > result_det(vec_size);
    std::vector< arma::mat::fixed<3,3> > result_mat33(vec_size);
    std::vector< arma::mat::fixed<3,2> > result_mat32(vec_size);

    START_TIMER("determinant_own");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_det[j] = determinant( mat33_vec[j] );
    END_TIMER("determinant_own");

    START_TIMER("determinant_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_det[j] = det( mat33_vec[j] );
    END_TIMER("determinant_arma");

    START_TIMER("inv_33_own");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat33[j] = inverse( mat33_vec[j] );
    END_TIMER("inv_33_own");

    START_TIMER("inv_33_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat33[j] = inv( mat33_vec[j] );
    END_TIMER("inv_33_arma");

    START_TIMER("pinv_23_own");
    for (uint i=0; i<N_RUNS/10; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat32[j] = inverse( mat23_vec[j] );
    END_TIMER("pinv_23_own");

    START_TIMER("pinv_23_arma");
    for (uint i=0; i<N_RUNS/10; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat32[j] = pinv( mat23_vec[j] );
    END_TIMER("pinv_23_arma");

    this->profiler_output("fem_tools");
}

