/*
 * fem_tools_bench.cpp
 *
 *  Created on: Oct 18, 2023
 *      Author: David Flanderka
 *
 *  Speed tests of determinant and inverse function
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
//#include <Eigen/QR>
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
TEST(FemToolsDevelopTest, functions_test) {
	using namespace fe_tools;

    arma::mat::fixed<1,1> mat11 = {2};
    arma::mat::fixed<2,2> mat22 = { {2, 3}, {4, 5} };
    arma::mat::fixed<3,3> mat33 = { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} };
    arma::mat::fixed<1,2> mat12 = { {2, 5} };
    arma::mat::fixed<2,1> mat21; mat21(0,0) = 2; mat21(1,0) = 5;
    arma::mat::fixed<1,3> mat13 = { {2, 4, 5} };
    arma::mat::fixed<3,1> mat31; mat31(0,0) = 2; mat31(1,0) = 4; mat31(2,0) = 5;
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
    arma::mat::fixed<2,1> inv12 = inverse(mat12);
    arma::mat::fixed<1,2> inv21 = inverse(mat21);
    arma::mat::fixed<3,1> inv13 = inverse(mat13);
    arma::mat::fixed<1,3> inv31 = inverse(mat31);
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
    // matrix 1x2
    arma::mat::fixed<1,1> multi_12 = mat12 * inv12;
    EXPECT_DOUBLE_EQ( multi_12(0,0), 1 );
    // matrix 3x1
    arma::mat::fixed<1,1> multi_21 = inv21 * mat21;
    EXPECT_DOUBLE_EQ( multi_21(0,0), 1 );
    // matrix 1x3
    arma::mat::fixed<1,1> multi_13 = mat13 * inv13;
    EXPECT_DOUBLE_EQ( multi_13(0,0), 1 );
    // matrix 3x1
    arma::mat::fixed<1,1> multi_31 = inv31 * mat31;
    EXPECT_DOUBLE_EQ( multi_31(0,0), 1 );
    // matrix 2x3
    arma::mat::fixed<2,2> multi_23 = mat23 * inv23;
    EXPECT_ARMA_EQ(expect_22, multi_23);
    // matrix 3x2
    arma::mat::fixed<2,2> multi_32 = inv32 * mat32;
    EXPECT_ARMA_EQ(expect_22, multi_32);
}


inline arma::vec::fixed<3> mat_multi_vec(const arma::mat::fixed<3,3> &mat, const arma::vec::fixed<3> &vec)
{
	arma::vec::fixed<3> res;

	res(0) = mat(0,0)*vec(0) + mat(0,1)*vec(1) + mat(0,2)*vec(2);
	res(1) = mat(1,0)*vec(0) + mat(1,1)*vec(1) + mat(1,2)*vec(2);
	res(2) = mat(2,0)*vec(0) + mat(2,1)*vec(1) + mat(2,2)*vec(2);

	return res;
}

//inline Eigen::Matrix<eigen_tools::Vec200,3,1> mat_multi_vec(const Eigen::Matrix<eigen_tools::Vec200,3,3> &mat, const Eigen::Matrix<eigen_tools::Vec200,3,1> &vec)
//{
//    Eigen::Matrix<eigen_tools::Vec200,3,1> res;
//
//    res(0) = mat(0,0)*vec(0) + mat(0,1)*vec(1) + mat(0,2)*vec(2);
//    res(1) = mat(1,0)*vec(0) + mat(1,1)*vec(1) + mat(1,2)*vec(2);
//    res(2) = mat(2,0)*vec(0) + mat(2,1)*vec(1) + mat(2,2)*vec(2);
//
//    return res;
//}



/**
 * Benchmark test. Compare speed of functions implemented in fem_tools.hh and armadillo library.
 *
 * Test compares following functions:
 *  - determinant (fem_tools and armadillo)
 *  - determinant (vectorized case in Armor object)
 *  - inversion of 3x3 matrix (fem_tools and armadillo)
 *  - pseudoinversion of 2x3 matrix (fem_tools and armadillo)
 *
 *
 *  Results:
 *   - result date: October 27, 2023:
 *   - run on: Dell Inspiron CPU 1.80 GHz, 16.0 GB RAM
 *   - n_repeats: 4e7
 *   - time unit: [s]
 *
 *               fem_tools   armadillo
 *  det 3x3         0.0899      0.4040
 *  inv 3x3         0.5873      1.5258
 *  pinv 2x3        1.9312     66.5879
 */
TEST_F(FemToolsTest, speed_test) {
	using namespace fe_tools;

    static const uint N_RUNS = 2e5;

    // Create arma and Eigen mats and vectors
    std::vector< arma::mat::fixed<3,3> > mat33_vec(200);      // vector of armadillo objects
    std::vector< arma::mat::fixed<2,3> > mat23_vec(200);
    std::vector< arma::vec::fixed<3> > vec3_vec(200);
//    Eigen::Matrix<eigen_tools::Vec200,3,1> vec3_eigen;             // Eigen matrix of vector items
//    Eigen::Matrix<eigen_tools::Vec200,2,3> mat23_eigen;
//    Eigen::Matrix<eigen_tools::Vec200,3,3> mat33_eigen;
//    std::vector< Eigen::Matrix<double,3,3> > eigen33_vec(200);  // vector of Eigen objects
//    std::vector< Eigen::Matrix<double,2,3> > eigen23_vec(200);
//    std::vector< Eigen::Matrix<double,3,1> > eigen3_vec(200);
    {
        // Fill arma and Eigen mats and vectors
        std::vector< arma::mat::fixed<3,3> > mat33_tmp = {
                { {1, 2, 3}, {2, 4, 5}, {3, 5, 6} },
                { {2, 4, 5}, {3, 2, 1}, {0, 6, 4} },
                { {1, 4, 2}, {5, 1, 3}, {2, 3, 4} },
                { {9, 7, 5}, {2, 4, 6}, {1, 3, 8} }
        };
        std::vector< arma::mat::fixed<2,3> > mat23_tmp = {
                { {1, 2, 3}, {2, 4, 5} },
                { {2, 4, 5}, {3, 2, 1} },
                { {1, 4, 2}, {5, 1, 3} },
                { {9, 7, 5}, {2, 4, 6} }
        };
        std::vector< arma::vec::fixed<3> > vec3_tmp = {
                {1, 2, 3},
                {4, 5, 6},
                {3, 2, 1},
                {6, 5, 4}
        };
        for (uint i=0; i<200; ++i) {
            uint data_i = i%4;
            mat33_vec[i] = mat33_tmp[data_i];
            mat23_vec[i] = mat23_tmp[data_i];
            vec3_vec[i] = vec3_tmp[data_i];
//            mat33_eigen(0,0)(i) = mat33_tmp[data_i](0,0);
//            mat33_eigen(0,1)(i) = mat33_tmp[data_i](0,1);
//            mat33_eigen(0,2)(i) = mat33_tmp[data_i](0,2);
//            mat33_eigen(1,0)(i) = mat33_tmp[data_i](1,0);
//            mat33_eigen(1,1)(i) = mat33_tmp[data_i](1,1);
//            mat33_eigen(1,2)(i) = mat33_tmp[data_i](1,2);
//            mat33_eigen(2,0)(i) = mat33_tmp[data_i](2,0);
//            mat33_eigen(2,1)(i) = mat33_tmp[data_i](2,1);
//            mat33_eigen(2,2)(i) = mat33_tmp[data_i](2,2);
//            mat23_eigen(0,0)(i) = mat23_tmp[data_i](0,0);
//            mat23_eigen(0,1)(i) = mat23_tmp[data_i](0,1);
//            mat23_eigen(0,2)(i) = mat23_tmp[data_i](0,2);
//            mat23_eigen(1,0)(i) = mat23_tmp[data_i](1,0);
//            mat23_eigen(1,1)(i) = mat23_tmp[data_i](1,1);
//            mat23_eigen(1,2)(i) = mat23_tmp[data_i](1,2);
//            vec3_eigen(0)(i) = vec3_tmp[data_i](0);
//            vec3_eigen(1)(i) = vec3_tmp[data_i](1);
//            vec3_eigen(2)(i) = vec3_tmp[data_i](2);
//            eigen33_vec[i](0,0) = mat33_tmp[data_i](0,0);
//            eigen33_vec[i](0,1) = mat33_tmp[data_i](0,1);
//            eigen33_vec[i](0,2) = mat33_tmp[data_i](0,2);
//            eigen33_vec[i](1,0) = mat33_tmp[data_i](1,0);
//            eigen33_vec[i](1,1) = mat33_tmp[data_i](1,1);
//            eigen33_vec[i](1,2) = mat33_tmp[data_i](1,2);
//            eigen33_vec[i](2,0) = mat33_tmp[data_i](2,0);
//            eigen33_vec[i](2,1) = mat33_tmp[data_i](2,1);
//            eigen33_vec[i](2,2) = mat33_tmp[data_i](2,2);
//            eigen23_vec[i](0,0) = mat23_tmp[data_i](0,0);
//            eigen23_vec[i](0,1) = mat23_tmp[data_i](0,1);
//            eigen23_vec[i](0,2) = mat23_tmp[data_i](0,2);
//            eigen23_vec[i](1,0) = mat23_tmp[data_i](1,0);
//            eigen23_vec[i](1,1) = mat23_tmp[data_i](1,1);
//            eigen23_vec[i](1,2) = mat23_tmp[data_i](1,2);
//            eigen3_vec[i](0) = vec3_tmp[data_i](0);
//            eigen3_vec[i](1) = vec3_tmp[data_i](1);
//            eigen3_vec[i](2) = vec3_tmp[data_i](2);
        }
    }

    // Declaration of result variables
    uint vec_size = mat33_vec.size();
    std::vector< double > result_det(vec_size);
    std::vector< arma::mat::fixed<3,3> > result_mat33(vec_size);
    std::vector< arma::mat::fixed<3,2> > result_mat32(vec_size);
    std::vector< arma::vec::fixed<3> > result_vec3(vec_size);
//    std::vector< double > eigen_result_det(vec_size);
//    std::vector< Eigen::Matrix<double,3,2> > eigen_result_mat32(vec_size);
//    std::vector< Eigen::Matrix<double,3,3> > eigen_result_mat33(vec_size);
//    std::vector< Eigen::Matrix<double,3,1> > eigen_result_vec3(vec_size);
//    eigen_tools::Vec200 result_det_eigen;
//    Eigen::Matrix<eigen_tools::Vec200,3,2> result_mat32_eigen;
//    Eigen::Matrix<eigen_tools::Vec200,3,3> result_mat33_eigen;
//    Eigen::Matrix<eigen_tools::Vec200,3,1> result_vec3_eigen;

    /**
     * Determinant 3x3
     */
    START_TIMER("DET");
    START_TIMER("determinant_fe_tools");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_det[j] = fe_tools::determinant( mat33_vec[j] );
    END_TIMER("determinant_fe_tools");

//    START_TIMER("determinant_eigen_tools");
//    for (uint i=0; i<N_RUNS; ++i)
//        result_det_eigen = eigen_tools::determinant( mat33_eigen );
//    END_TIMER("determinant_eigen_tools");
//
//    START_TIMER("determinant_eigen_vec");
//    for (uint i=0; i<N_RUNS; ++i)
//        result_det_eigen = mat33_eigen.determinant();
//    END_TIMER("determinant_eigen_vec");
//
//    START_TIMER("determinant_eigen");
//    for (uint i=0; i<N_RUNS; ++i)
//        for (uint j=0; j<vec_size; ++j) eigen_result_det[j] = eigen33_vec[j].determinant();
//    END_TIMER("determinant_eigen");

    START_TIMER("determinant_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_det[j] = det( mat33_vec[j] );
    END_TIMER("determinant_arma");
    END_TIMER("DET");

    /**
     * Inverse 3x3
     */
    START_TIMER("INV");
    START_TIMER("inv_33_fe_tools");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat33[j] = fe_tools::inverse( mat33_vec[j] );
    END_TIMER("inv_33_fe_tools");

//    START_TIMER("inv_33_eigen_tools");
//    for (uint i=0; i<N_RUNS; ++i)
//        result_mat33_eigen = eigen_tools::inverse<3,3>( mat33_eigen );
//    END_TIMER("inv_33_eigen_tools");
//
//    START_TIMER("inv_33_eigen_vec");
//    for (uint i=0; i<N_RUNS; ++i)
//        result_mat33_eigen = mat33_eigen.inverse();
//    END_TIMER("inv_33_eigen_vec");
//
//    START_TIMER("inv_33_eigen");
//    for (uint i=0; i<N_RUNS; ++i)
//        for (uint j=0; j<vec_size; ++j) eigen_result_mat33[j] = eigen33_vec[j].inverse();
//    END_TIMER("inv_33_eigen");

    START_TIMER("inv_33_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat33[j] = inv( mat33_vec[j] );
    END_TIMER("inv_33_arma");
    END_TIMER("INV");

    /**
     * Pseudoinverse 2x3
     */
    START_TIMER("PINV");
    START_TIMER("pinv_23_fe_tools");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat32[j] = fe_tools::inverse( mat23_vec[j] );
    END_TIMER("pinv_23_fe_tools");

//    START_TIMER("pinv_23_eigen_tools");
//    for (uint i=0; i<N_RUNS; ++i)
//        result_mat32_eigen = eigen_tools::inverse<2,3>( mat23_eigen );
//    END_TIMER("pinv_23_eigen_tools");
//
//    START_TIMER("pinv_23_eigen_vec");
//    // https://stackoverflow.com/questions/44465197/eigen-library-pseudo-inverse-of-matrix-matlab-pinv
//    for (uint i=0; i<N_RUNS; ++i)
//        result_mat32_eigen = mat23_eigen.completeOrthogonalDecomposition().pseudoInverse();
//    END_TIMER("pinv_23_eigen_vec");
//
//    START_TIMER("pinv_23_eigen");
//    for (uint i=0; i<N_RUNS; ++i)
//        for (uint j=0; j<vec_size; ++j) eigen_result_mat32[j] = eigen23_vec[j].completeOrthogonalDecomposition().pseudoInverse();
//    END_TIMER("pinv_23_eigen");

    START_TIMER("pinv_23_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_mat32[j] = pinv( mat23_vec[j] );
    END_TIMER("pinv_23_arma");
    END_TIMER("PINV");

    /**
     * mat3x3 @ vec3
     */
    START_TIMER("MULTI");
    START_TIMER("multi_mat_vec_fe_tools");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_vec3[j] = mat_multi_vec(mat33_vec[j], vec3_vec[j]);
    END_TIMER("multi_mat_vec_fe_tools");

//    START_TIMER("multi_mat_vec_eigen_tools");
//    for (unsigned int i=0; i<N_RUNS; ++i) {
//    	result_vec3_eigen = mat_multi_vec(mat33_eigen, vec3_eigen);
//    }
//    END_TIMER("multi_mat_vec_eigen_tools");
//
//    START_TIMER("multi_mat_vec_eigen_vec");
//    for (unsigned int i=0; i<N_RUNS; ++i) {
//    	result_vec3_eigen = mat33_eigen * vec3_eigen;
//    }
//    END_TIMER("multi_mat_vec_eigen_vec");
//
//    START_TIMER("multi_mat_vec_eigen");
//    for (uint i=0; i<N_RUNS; ++i)
//        for (uint j=0; j<vec_size; ++j) eigen_result_vec3[j] = eigen33_vec[j] * eigen3_vec[j];
//    END_TIMER("multi_mat_vec_eigen");

    START_TIMER("multi_mat_vec_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_vec3[j] = mat33_vec[j] * vec3_vec[j];
    END_TIMER("multi_mat_vec_arma");
    END_TIMER("MULTI");

    this->profiler_output("fem_tools");
}

