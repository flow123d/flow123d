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
#include <armadillo>
#include "arma_expect.hh"
#include <Eigen/QR>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "fem/fem_tools.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"
#include "system/armor.hh"
#include "system/fmt/posix.h"           // for FMT_UNUSED


//template <unsigned int size>
//class VectorCol {
//public:
//    typename Eigen::Array<double,size,1> data;
//
//    VectorCol() {}
//
//    VectorCol(FMT_UNUSED int n) {}
//
//    inline double & operator[](std::size_t item) {
//        return data[item];
//    }
//
//    inline double & operator()(std::size_t item) {
//        return data[item];
//    }
//
//    inline VectorCol<size> operator+(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//        res.data = this->data + other.data;
//        return res;
//    }
//
//    /// Binary operator minus
//    inline VectorCol<size> operator-(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//        res.data = this->data - other.data;
//        return res;
//    }
//
//    /// Unary operator minus
//    inline VectorCol<size> operator-() const {
//    	VectorCol<size> res;
//        res.data = - this->data;
//        return res;
//    }
//
//    inline VectorCol<size> operator*(const double &coef) const {
//    	VectorCol<size> res;
//        res.data = this->data * coef;
//        return res;
//    }
//
//    inline VectorCol<size> operator/(const double &coef) const {
//    	VectorCol<size> res;
//        res.data = this->data / coef;
//        return res;
//    }
//
//    inline VectorCol<size> operator*(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//    	res.data = this->data * other.data;
//        return res;
//    }
//
//    inline VectorCol<size> operator/(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//    	res.data = this->data / other.data;
//        return res;
//    }
//
//    inline VectorCol<size> inverse() const {
//    	VectorCol<size> res;
//    	res.data = 1 / this->data;
//        return res;
//    }
//
//    inline VectorCol<size> sqrt() const {
//    	VectorCol<size> res;
//    	res.data = this->data.abs().sqrt();
//        return res;
//    }
//
//};
//
//template class VectorCol<8>;
//template class VectorCol<200>;
//
//
//namespace Eigen {
//    template<> struct NumTraits<VectorCol<8>> : GenericNumTraits<VectorCol<8>>
//    {
//        typedef VectorCol<8> Real;
//        typedef VectorCol<8> NonInteger;
//        typedef VectorCol<8> Nested;
//
//        enum {
//            IsInteger = 0,
//            IsSigned = 1,
//            IsComplex = 0,
//            RequireInitialization = 1,
//            ReadCost = 6,
//            AddCost = 150,
//            MulCost = 100
//        };
//    };
//    template<> struct NumTraits<VectorCol<200>> : GenericNumTraits<VectorCol<200>>
//    {
//        typedef VectorCol<200> Real;
//        typedef VectorCol<200> NonInteger;
//        typedef VectorCol<200> Nested;
//
//        enum {
//            IsInteger = 0,
//            IsSigned = 1,
//            IsComplex = 0,
//            RequireInitialization = 1,
//            ReadCost = 6,
//            AddCost = 150,
//            MulCost = 100
//        };
//    };
//}


namespace eigen_tools_test {

// use only one of the following typedef
//typedef VectorCol<200> UserVec;
typedef Eigen::Array<double,Eigen::Dynamic,1> UserVec;



/**
 * @brief Calculates determinant of a rectangular matrix.
 */
template<class T>
UserVec determinant(const T &M);



inline Eigen::Matrix<UserVec,1,1> normal_matrix(const Eigen::Matrix<UserVec,1,2> &A) {
	Eigen::Matrix<UserVec,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1);
    return res;
}

inline Eigen::Matrix<UserVec,1,1> normal_matrix(const Eigen::Matrix<UserVec,2,1> &A) {
	Eigen::Matrix<UserVec,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0);
    return res;
}

inline Eigen::Matrix<UserVec,1,1> normal_matrix(const Eigen::Matrix<UserVec,1,3> &A) {
	Eigen::Matrix<UserVec,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    return res;
}

inline Eigen::Matrix<UserVec,1,1> normal_matrix(const Eigen::Matrix<UserVec,3,1> &A) {
	Eigen::Matrix<UserVec,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    return res;
}

inline Eigen::Matrix<UserVec,2,2> normal_matrix(const Eigen::Matrix<UserVec,2,3> &A) {
    Eigen::Matrix<UserVec,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    res(0,1) = A(0,0)*A(1,0)+A(0,1)*A(1,1)+A(0,2)*A(1,2);
    res(1,0) = A(1,0)*A(0,0)+A(1,1)*A(0,1)+A(1,2)*A(0,2);
    res(1,1) = A(1,0)*A(1,0)+A(1,1)*A(1,1)+A(1,2)*A(1,2);
    return res;
}

inline Eigen::Matrix<UserVec,2,2> normal_matrix(const Eigen::Matrix<UserVec,3,2> &A) {
	Eigen::Matrix<UserVec,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    res(0,1) = A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    res(1,0) = A(0,1)*A(0,0)+A(1,1)*A(1,0)+A(2,1)*A(2,0);
    res(1,1) = A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    return res;
}



template<> inline UserVec determinant(const Eigen::Matrix<UserVec,1,1> &M)
{
    return M(0,0);
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,2,2> &M)
{
    return M(0,0)*M(1,1) - M(1,0)*M(0,1);
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,3,3> &M)
{
    return ( M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) )
         - ( M(2,0)*M(1,1)*M(0,2) + M(2,1)*M(1,2)*M(0,0) + M(2,2)*M(1,0)*M(0,1) );
}

template<> inline UserVec determinant(FMT_UNUSED const Eigen::Matrix<UserVec,0,3> &M)
{
    UserVec vec; vec.resize(200);
    return vec;
}

template<> inline UserVec determinant(FMT_UNUSED const Eigen::Matrix<UserVec,3,0> &M)
{
    UserVec vec; vec.resize(200);
    return vec;
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,1,2> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,2,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,1,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,3,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,2,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline UserVec determinant(const Eigen::Matrix<UserVec,3,2> &M)
{
	return determinant(normal_matrix(M)).sqrt();
}


/**
 * @brief Calculates inverse of rectangular matrix or pseudoinverse of non-rectangular matrix.
 */
template<int m, int n>
Eigen::Matrix<UserVec,n,m> inverse(const Eigen::Matrix<UserVec,m,n> &A) {
    // only for cases m > n
    return inverse(normal_matrix(A)) * A.transpose();
}


template<> inline Eigen::Matrix<UserVec,1,1> inverse<1,1>(const Eigen::Matrix<UserVec,1,1> &A)
{
	Eigen::Matrix<UserVec,1,1> B;
    B(0,0) = A(0,0).inverse(); // 1/A(0,0)
    return B;
}

template<> inline Eigen::Matrix<UserVec,2,2> inverse<2,2>(const Eigen::Matrix<UserVec,2,2> &A)
{
	Eigen::Matrix<UserVec,2,2> B;
	UserVec det = determinant(A);

    B(0,0) = A(1,1) / det;
    B(0,1) = -A(0,1) / det;
    B(1,0) = -A(1,0) / det;
    B(1,1) = A(0,0) / det;
    return B;
}

template<> inline Eigen::Matrix<UserVec,3,3> inverse<3,3>(const Eigen::Matrix<UserVec,3,3> &A)
{
    Eigen::Matrix<UserVec,3,3> B;

    B(0,0) = A(1,1)*A(2,2) - A(2,1)*A(1,2);
    B(0,1) = A(2,0)*A(1,2) - A(1,0)*A(2,2);
    B(0,2) = A(1,0)*A(2,1) - A(2,0)*A(1,1);

    UserVec det = A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
    B(0,0) = B(0,0) / det;
    B(0,1) = B(0,1) / det;
    B(0,2) = B(0,2) / det;

    B(1,0) = (A(2,1)*A(0,2) - A(0,1)*A(2,2)) / det;
    B(1,1) = (A(0,0)*A(2,2) - A(2,0)*A(0,2)) / det;
    B(1,2) = (A(2,0)*A(0,1) - A(0,0)*A(2,1)) / det;

    B(2,0) = (A(0,1)*A(1,2) - A(1,1)*A(0,2)) / det;
    B(2,1) = (A(1,0)*A(0,2) - A(0,0)*A(1,2)) / det;
    B(2,2) = (A(0,0)*A(1,1) - A(1,0)*A(0,1)) / det;

    return B;
}

template<> inline Eigen::Matrix<UserVec,2,1> inverse<1,2>(const Eigen::Matrix<UserVec,1,2> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<UserVec,3,1> inverse<1,3>(const Eigen::Matrix<UserVec,1,3> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<UserVec,3,2> inverse<2,3>(const Eigen::Matrix<UserVec,2,3> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}


} // closing namespace eigen_tools_test



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
    arma::mat::fixed<3,3> mat33 = { {1, 2, 3}, {5, 4, 5}, {2, 3, 6} };
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

inline Eigen::Matrix<eigen_tools_test::UserVec,3,1> mat_multi_vec(const Eigen::Matrix<eigen_tools_test::UserVec,3,3> &mat, const Eigen::Matrix<eigen_tools_test::UserVec,3,1> &vec)
{
    Eigen::Matrix<eigen_tools_test::UserVec,3,1> res;

    res(0) = mat(0,0)*vec(0) + mat(0,1)*vec(1) + mat(0,2)*vec(2);
    res(1) = mat(1,0)*vec(0) + mat(1,1)*vec(1) + mat(1,2)*vec(2);
    res(2) = mat(2,0)*vec(0) + mat(2,1)*vec(1) + mat(2,2)*vec(2);

    return res;
}



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
    std::vector< Eigen::Matrix<double,3,3> > eigen33_vec(200);  // vector of Eigen objects
    std::vector< Eigen::Matrix<double,2,3> > eigen23_vec(200);
    std::vector< Eigen::Matrix<double,3,1> > eigen3_vec(200);
    Eigen::Matrix<eigen_tools_test::UserVec,3,1> vec3_eigen;             // Eigen matrix of vector items
    Eigen::Matrix<eigen_tools_test::UserVec,2,3> mat23_eigen;
    Eigen::Matrix<eigen_tools_test::UserVec,3,3> mat33_eigen;
    for (uint i=0; i<vec3_eigen.size(); ++i)
        vec3_eigen(i).resize(200);
    for (uint r=0; r<mat23_eigen.rows(); ++r)
        for (uint c=0; c<mat23_eigen.cols(); ++c)
            mat23_eigen(r,c).resize(200);
    for (uint r=0; r<mat33_eigen.rows(); ++r)
        for (uint c=0; c<mat33_eigen.cols(); ++c)
            mat33_eigen(r,c).resize(200);
    {
        // Fill arma and Eigen mats and vectors
        std::vector< arma::mat::fixed<3,3> > mat33_tmp = {
                { {1, 2, 3}, {5, 4, 5}, {2, 3, 6} },
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
            mat33_eigen(0,0)(i) = mat33_tmp[data_i](0,0);
            mat33_eigen(0,1)(i) = mat33_tmp[data_i](0,1);
            mat33_eigen(0,2)(i) = mat33_tmp[data_i](0,2);
            mat33_eigen(1,0)(i) = mat33_tmp[data_i](1,0);
            mat33_eigen(1,1)(i) = mat33_tmp[data_i](1,1);
            mat33_eigen(1,2)(i) = mat33_tmp[data_i](1,2);
            mat33_eigen(2,0)(i) = mat33_tmp[data_i](2,0);
            mat33_eigen(2,1)(i) = mat33_tmp[data_i](2,1);
            mat33_eigen(2,2)(i) = mat33_tmp[data_i](2,2);
            mat23_eigen(0,0)(i) = mat23_tmp[data_i](0,0);
            mat23_eigen(0,1)(i) = mat23_tmp[data_i](0,1);
            mat23_eigen(0,2)(i) = mat23_tmp[data_i](0,2);
            mat23_eigen(1,0)(i) = mat23_tmp[data_i](1,0);
            mat23_eigen(1,1)(i) = mat23_tmp[data_i](1,1);
            mat23_eigen(1,2)(i) = mat23_tmp[data_i](1,2);
            vec3_eigen(0)(i) = vec3_tmp[data_i](0);
            vec3_eigen(1)(i) = vec3_tmp[data_i](1);
            vec3_eigen(2)(i) = vec3_tmp[data_i](2);
            eigen33_vec[i](0,0) = mat33_tmp[data_i](0,0);
            eigen33_vec[i](0,1) = mat33_tmp[data_i](0,1);
            eigen33_vec[i](0,2) = mat33_tmp[data_i](0,2);
            eigen33_vec[i](1,0) = mat33_tmp[data_i](1,0);
            eigen33_vec[i](1,1) = mat33_tmp[data_i](1,1);
            eigen33_vec[i](1,2) = mat33_tmp[data_i](1,2);
            eigen33_vec[i](2,0) = mat33_tmp[data_i](2,0);
            eigen33_vec[i](2,1) = mat33_tmp[data_i](2,1);
            eigen33_vec[i](2,2) = mat33_tmp[data_i](2,2);
            eigen23_vec[i](0,0) = mat23_tmp[data_i](0,0);
            eigen23_vec[i](0,1) = mat23_tmp[data_i](0,1);
            eigen23_vec[i](0,2) = mat23_tmp[data_i](0,2);
            eigen23_vec[i](1,0) = mat23_tmp[data_i](1,0);
            eigen23_vec[i](1,1) = mat23_tmp[data_i](1,1);
            eigen23_vec[i](1,2) = mat23_tmp[data_i](1,2);
            eigen3_vec[i](0) = vec3_tmp[data_i](0);
            eigen3_vec[i](1) = vec3_tmp[data_i](1);
            eigen3_vec[i](2) = vec3_tmp[data_i](2);
        }
    }

    // Declaration of result variables
    uint vec_size = mat33_vec.size();
    std::vector< double > result_det(vec_size);
    std::vector< arma::mat::fixed<3,3> > result_mat33(vec_size);
    std::vector< arma::mat::fixed<3,2> > result_mat32(vec_size);
    std::vector< arma::vec::fixed<3> > result_vec3(vec_size);
    std::vector< double > eigen_result_det(vec_size);
    std::vector< Eigen::Matrix<double,3,2> > eigen_result_mat32(vec_size);
    std::vector< Eigen::Matrix<double,3,3> > eigen_result_mat33(vec_size);
    std::vector< Eigen::Matrix<double,3,1> > eigen_result_vec3(vec_size);
    eigen_tools_test::UserVec result_det_eigen;
    Eigen::Matrix<eigen_tools_test::UserVec,3,2> result_mat32_eigen;
    Eigen::Matrix<eigen_tools_test::UserVec,3,3> result_mat33_eigen;
    Eigen::Matrix<eigen_tools_test::UserVec,3,1> result_vec3_eigen;

    /**
     * Determinant 3x3
     */
    START_TIMER("DET");
    START_TIMER("determinant_fe_tools");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_det[j] = fe_tools::determinant( mat33_vec[j] );
    END_TIMER("determinant_fe_tools");

    START_TIMER("determinant_eigen_tools");
    for (uint i=0; i<N_RUNS; ++i)
        result_det_eigen = eigen_tools_test::determinant( mat33_eigen );
    END_TIMER("determinant_eigen_tools");

    START_TIMER("determinant_eigen_vec");
    for (uint i=0; i<N_RUNS; ++i)
        result_det_eigen = mat33_eigen.determinant();
    END_TIMER("determinant_eigen_vec");

    START_TIMER("determinant_eigen");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) eigen_result_det[j] = eigen33_vec[j].determinant();
    END_TIMER("determinant_eigen");

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

    START_TIMER("inv_33_eigen_tools");
    for (uint i=0; i<N_RUNS; ++i)
        result_mat33_eigen = eigen_tools_test::inverse<3,3>( mat33_eigen );
    END_TIMER("inv_33_eigen_tools");

    START_TIMER("inv_33_eigen_vec");
    for (uint i=0; i<N_RUNS; ++i)
        result_mat33_eigen = mat33_eigen.inverse();
    END_TIMER("inv_33_eigen_vec");

    START_TIMER("inv_33_eigen");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) eigen_result_mat33[j] = eigen33_vec[j].inverse();
    END_TIMER("inv_33_eigen");

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

    START_TIMER("pinv_23_eigen_tools");
    for (uint i=0; i<N_RUNS; ++i)
        result_mat32_eigen = eigen_tools_test::inverse<2,3>( mat23_eigen );
    END_TIMER("pinv_23_eigen_tools");

//    START_TIMER("pinv_23_eigen_vec");
//    // https://stackoverflow.com/questions/44465197/eigen-library-pseudo-inverse-of-matrix-matlab-pinv
//    for (uint i=0; i<N_RUNS; ++i)
//        result_mat32_eigen = mat23_eigen.completeOrthogonalDecomposition().pseudoInverse();
//    END_TIMER("pinv_23_eigen_vec");

    START_TIMER("pinv_23_eigen");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) eigen_result_mat32[j] = eigen23_vec[j].completeOrthogonalDecomposition().pseudoInverse();
    END_TIMER("pinv_23_eigen");

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

    START_TIMER("multi_mat_vec_eigen_tools");
    for (unsigned int i=0; i<N_RUNS; ++i) {
    	result_vec3_eigen = mat_multi_vec(mat33_eigen, vec3_eigen);
    }
    END_TIMER("multi_mat_vec_eigen_tools");

    START_TIMER("multi_mat_vec_eigen_vec");
    for (unsigned int i=0; i<N_RUNS; ++i) {
    	result_vec3_eigen = mat33_eigen * vec3_eigen;
    }
    END_TIMER("multi_mat_vec_eigen_vec");

    START_TIMER("multi_mat_vec_eigen");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) eigen_result_vec3[j] = eigen33_vec[j] * eigen3_vec[j];
    END_TIMER("multi_mat_vec_eigen");

    START_TIMER("multi_mat_vec_arma");
    for (uint i=0; i<N_RUNS; ++i)
        for (uint j=0; j<vec_size; ++j) result_vec3[j] = mat33_vec[j] * vec3_vec[j];
    END_TIMER("multi_mat_vec_arma");
    END_TIMER("MULTI");

    this->profiler_output("fem_tools");
}

#include "fem/eigen_tools.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"

/// Compute used space and free space of arena in operations
//TEST_F(FemToolsTest, free_space_test) {
//    static const uint data_size = 30;
//
//    AssemblyArena asm_arena(1024 * 1024, 256);
//
//    // Defines Eigen objects
//    // Speed test - comparation ArrayDbl and ArenaVec - inputs
//    Eigen::Matrix<ArrayDbl,3,3> mat33_array;
//    Eigen::Matrix<ArrayDbl,2,3> mat23_array;
//    Eigen::Matrix<ArrayDbl,3,1> vec3_array;
//    Eigen::Matrix<ArenaVec<double>, 3, 3> mat33_arena;
//    Eigen::Matrix<ArenaVec<double>, 2, 3> mat23_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 1> vec3_arena;
//    // Speed test - comparation ArrayDbl and ArenaVec - outputs
//    ArrayDbl result_det_array;
//    Eigen::Matrix<ArrayDbl,3,3> result_inv_array;
//    Eigen::Matrix<ArrayDbl,3,2> result_pinv_array;
//    Eigen::Matrix<ArrayDbl,3,1> result_multi_array;
//    ArenaVec<double> result_det_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 3> result_inv_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 2> result_pinv_arena;
//    Eigen::Matrix<ArenaVec<double>, 3, 1> result_multi_arena;
//
//    // Fills Eigen objects
//    for (uint i=0; i<9; ++i) {
//    	mat33_arena(i) = ArenaVec<double>(data_size, asm_arena);
//        for (uint j=0; j<data_size; ++j)
//            mat33_arena(i)(j) = 0.1 * (2*i + j + 1);
//    }
//    for (uint i=0; i<6; ++i) {
//    	mat23_arena(i) = ArenaVec<double>(data_size, asm_arena);
//        for (uint j=0; j<data_size; ++j)
//            mat23_arena(i)(j) = 0.1 * (2*i + j + 1);
//    }
//    for (uint i=0; i<3; ++i) {
//    	vec3_arena(i) = ArenaVec<double>(data_size, asm_arena);
//        for (uint j=0; j<data_size; ++j)
//            vec3_arena(i)(j) = 0.1 * (2*i + j + 1);
//    }
//
//    result_det_arena = eigen_arena_tools::determinant<3,3>(mat33_arena);
//    result_inv_arena = eigen_arena_tools::inverse<3,3>(mat33_arena);
//    result_pinv_arena = eigen_arena_tools::inverse<2,3>(mat23_arena);
//    result_multi_arena = mat33_arena * vec3_arena;
//    asm_arena.print_space();
//
//}


//class OpDeterminant {
//public>
//    inline ArrayDbl test_array_dbl(Eigen::Matrix<ArrayDbl,3,3> &mat33_array) {
//        return eigen_tools::determinant<3,3>(mat33_array);
//    }
//
//    inline ArenaVec<double> test_arena(Eigen::Matrix<ArenaVec<double>, 3, 3> &mat33_arena) {
//        return eigen_arena_tools::determinant<3,3>(mat33_arena);
//    }
//};



class ArenaSpeedTest {
public:
    static const uint N_RUNS = 4e5;

    /// Constructor, pass arena size in kB
    ArenaSpeedTest(size_t arena_size, size_t simd_size = 32)
    : asm_arena(arena_size * 1024, simd_size)
    {}

    /// Fill mat33_array, mat33_arena
    void fill_matt33(uint data_size) {
        for (uint i=0; i<9; ++i) {
            mat33_array(i).resize(data_size);
            mat33_arena(i) = ArenaVec<double>(data_size, asm_arena);
            for (uint j=0; j<data_size; ++j) {
                mat33_array(i)(j) = 0.1 * (2*i + j + 1);
                mat33_arena(i)(j) = 0.1 * (2*i + j + 1);
            }
        }
    }

    /// Same for matrix 2x3
    void fill_matt23(uint data_size) {
        for (uint i=0; i<6; ++i) {
            mat23_array(i).resize(data_size);
            mat23_arena(i) = ArenaVec<double>(data_size, asm_arena);
            for (uint j=0; j<data_size; ++j) {
                mat23_array(i)(j) = 0.1 * (2*i + j + 1);
                mat23_arena(i)(j) = 0.1 * (2*i + j + 1);
            }
        }
    }

    /// Same for vector 3
    void fill_vec3(uint data_size) {
        for (uint i=0; i<3; ++i) {
            vec3_array(i).resize(data_size);
            vec3_arena(i) = ArenaVec<double>(data_size, asm_arena);
            for (uint j=0; j<data_size; ++j) {
                vec3_array(i)(j) = 0.1 * (2*i + j + 1);
                vec3_arena(i)(j) = 0.1 * (2*i + j + 1);
            }
        }
    }

    /// Set patch arena to mat33_arena
    void set_patch_arena_mat33() {
        for (uint i=0; i<9; ++i)
            mat33_arena(i).set_patch_arena(*patch_arena);
    }

    /// Same for matrix 2x3
    void set_patch_arena_mat23() {
        for (uint i=0; i<6; ++i)
            mat23_arena(i).set_patch_arena(*patch_arena);
    }

    /// Same for vector 3
    void set_patch_arena_vec3() {
        for (uint i=0; i<3; ++i)
            vec3_arena(i).set_patch_arena(*patch_arena);
    }

    /**
     * Determinant Matrix3x3
     */
    void determinant_test(uint data_size) {
        fill_matt33(data_size);
        patch_arena = asm_arena.get_child_arena();
        set_patch_arena_mat33();

        START_TIMER("determinant_array_dbl");
        for (uint i=0; i<N_RUNS; ++i)
            result_det_array = eigen_tools::determinant<3,3>(mat33_array);
        END_TIMER("determinant_array_dbl");

        START_TIMER("determinant_arena_vec");
        for (uint i=0; i<N_RUNS; ++i) {
            patch_arena->reset();
            result_det_arena = eigen_arena_tools::determinant<3,3>(mat33_arena);
        }
        END_TIMER("determinant_arena_vec");
        EXPECT_DOUBLE_EQ( result_det_array(0), result_det_arena(0) );
    }

    /**
     * Inverse Matrix3x3
     */
    void inverse3x3_test(uint data_size) {
        fill_matt33(data_size);
        patch_arena = asm_arena.get_child_arena();
        set_patch_arena_mat33();

        START_TIMER("inverse_array_dbl");
        for (uint i=0; i<N_RUNS; ++i)
            result_inv_array = eigen_tools::inverse<3,3>(mat33_array);
        END_TIMER("inverse_array_dbl");

        START_TIMER("inverse_arena_vec");
        for (uint i=0; i<N_RUNS; ++i) {
            patch_arena->reset();
            result_inv_arena = eigen_arena_tools::inverse<3,3>(mat33_arena);
        }
        END_TIMER("inverse_arena_vec");
        EXPECT_DOUBLE_EQ( result_inv_array(0)(0), result_inv_arena(0)(0) );
    }

    /**
     * Pseudoinverse Matrix2x3
     */
    void inverse2x3_test(uint data_size) {
        fill_matt23(data_size);
        patch_arena = asm_arena.get_child_arena();
        set_patch_arena_mat23();

        START_TIMER("pinverse_array_dbl");
        for (uint i=0; i<N_RUNS; ++i)
            result_pinv_array = eigen_tools::inverse<2,3>(mat23_array);
        END_TIMER("pinverse_array_dbl");

        START_TIMER("pinverse_arena_vec");
        for (uint i=0; i<N_RUNS; ++i) {
            patch_arena->reset();
            result_pinv_arena = eigen_arena_tools::inverse<2,3>(mat23_arena);
        }
        END_TIMER("pinverse_arena_vec");
        EXPECT_DOUBLE_EQ( result_pinv_array(0)(0), result_pinv_arena(0)(0) );
    }

    /**
     * Multiplication Matrix3x3 * Vector3
     */
    void mat_multi_test(uint data_size) {
    	fill_matt33(data_size);
        fill_vec3(data_size);
        patch_arena = asm_arena.get_child_arena();
        set_patch_arena_mat33();
        set_patch_arena_vec3();

        START_TIMER("multi_array_dbl");
        for (uint i=0; i<N_RUNS; ++i)
        	result_multi_array = mat33_array * vec3_array;
        END_TIMER("multi_array_dbl");

        START_TIMER("multi_arena_vec");
        for (uint i=0; i<N_RUNS; ++i) {
            patch_arena->reset();
            result_multi_arena = mat33_arena * vec3_arena;
        }
        END_TIMER("multi_arena_vec");
        EXPECT_DOUBLE_EQ( result_multi_array(0)(0), result_multi_arena(0)(0) );
    }

    /**
     * Test of all functions on given arena
     */
    void all_tests(std::vector<uint> data_sizes) {
        determinant_test(data_sizes[0]);
        reset_arena();
        inverse3x3_test(data_sizes[1]);
        reset_arena();
        inverse2x3_test(data_sizes[2]);
        reset_arena();
        mat_multi_test(data_sizes[3]);
        reset_arena();
    }

    /// Reset asm_arena and destroy patch_arena if exists.
    void reset_arena() {
    	if (patch_arena != nullptr) {
    	    delete patch_arena;
    	    patch_arena = nullptr;
    	}
    	asm_arena.reset();
    }

    // arena objects
    AssemblyArena asm_arena;
    PatchArena *patch_arena;

    // inputs
    Eigen::Matrix<ArrayDbl,3,3> mat33_array;
    Eigen::Matrix<ArrayDbl,2,3> mat23_array;
    Eigen::Matrix<ArrayDbl,3,1> vec3_array;
    Eigen::Matrix<ArenaVec<double>, 3, 3> mat33_arena;
    Eigen::Matrix<ArenaVec<double>, 2, 3> mat23_arena;
    Eigen::Matrix<ArenaVec<double>, 3, 1> vec3_arena;
    // results
    ArrayDbl result_det_array;
    Eigen::Matrix<ArrayDbl,3,3> result_inv_array;
    Eigen::Matrix<ArrayDbl,3,2> result_pinv_array;
    Eigen::Matrix<ArrayDbl,3,1> result_multi_array;
    ArenaVec<double> result_det_arena;
    Eigen::Matrix<ArenaVec<double>, 3, 3> result_inv_arena;
    Eigen::Matrix<ArenaVec<double>, 3, 2> result_pinv_arena;
    Eigen::Matrix<ArenaVec<double>, 3, 1> result_multi_arena;
};

void all_tests(size_t arena_size, std::vector<uint> data_sizes) {
    ArenaSpeedTest test_det(arena_size);
    test_det.determinant_test(data_sizes[0]);

    ArenaSpeedTest test_inv(arena_size);
    test_inv.inverse3x3_test(data_sizes[1]);

    ArenaSpeedTest test_pinv(arena_size);
    test_pinv.inverse2x3_test(data_sizes[2]);

    ArenaSpeedTest test_multi(arena_size);
    test_multi.mat_multi_test(data_sizes[3]);
}


TEST_F(FemToolsTest, speed_eigen_test) {
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> OuterDbl;

    MessageOut() << "L1 half ... " << std::endl;
    { // L1 half
        START_TIMER("L1_HALF");
        all_tests( 32, {156, 80, 76, 148} );
        END_TIMER("L1_HALF");
    }
    MessageOut() << " ...  OK" << std::endl;

    MessageOut() << "L1 full ... " << std::endl;
    { // L1 full
        START_TIMER("L1_FULL");
        all_tests( 64, {312, 160, 152, 300} );
        END_TIMER("L1_FULL");
    }
    MessageOut() << " ...  OK" << std::endl;

    MessageOut() << "L2 half ... " << std::endl;
    { // L2 half
        START_TIMER("L2_HALF");
        all_tests( 128, {628, 324, 308, 604} );
        END_TIMER("L2_HALF");
    }
    MessageOut() << " ...  OK" << std::endl;

    MessageOut() << "L2 full ... " << std::endl;
    { // L2 full
        START_TIMER("L2_FULL");
        all_tests(256, {1260, 652, 616, 1212} );
        END_TIMER("L2_FULL");
    }
    MessageOut() << " ...  OK" << std::endl;

    MessageOut() << "L3 half ... " << std::endl;
    { // L3 half
        START_TIMER("L3_HALF");
        all_tests(4 * 1024, {20164, 10484, 9892, 19416} );
        END_TIMER("L3_HALF");
    }
    MessageOut() << " ...  OK" << std::endl;

    MessageOut() << "L3 full ... " << std::endl;
    { // L3 full
        START_TIMER("L3_FULL");
        all_tests(8 * 1024, {40328, 20968, 19784, 38836} );
        END_TIMER("L3_FULL");
    }
    MessageOut() << " ...  OK" << std::endl;

    /* OUTER PRODUCT - compare speed of transpose and map */
    static const uint N_RUNS = 4e5;
    static const uint data_size = 30;

    // Speed test of outer product
    OuterDbl outer_prod1;
    OuterDbl outer_prod2;
    OuterDbl result_oprod_transp;
    OuterDbl result_oprod_map;

    // Fills Eigen objects
    outer_prod1.resize(3);
    for (uint i=0; i<3; ++i)
        outer_prod1(i) = 1 + i*0.1;
    outer_prod2.resize(data_size/3);
    for (uint i=0; i<data_size/3; ++i)
        outer_prod2(i) = 2.1 + i*0.2;
    result_oprod_transp.resize(data_size);
    result_oprod_map.resize(data_size);

    /**
     * OuterProduct compare speed of:
     *    ColVec3 * ColVec10.transpose
     *    ColVec3 * RowVec10 (RowVec10 is created as mapped ColVec10)
     */
    START_TIMER("OPROD");
    START_TIMER("outer_prod_transpose");
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> result_transp_ref(result_oprod_transp.data(), 3, data_size/3);
    for (uint i=0; i<N_RUNS; ++i)
    	result_transp_ref = outer_prod1 * outer_prod2.transpose();
    END_TIMER("outer_prod_transpose");

    START_TIMER("outer_prod_map");
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> result_map_ref(result_oprod_map.data(), 3, data_size/3);
    for (uint i=0; i<N_RUNS; ++i) {
        Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic>> outer_prod2_ref(outer_prod2.data(), 1, data_size/3);
        result_map_ref = outer_prod1 * outer_prod2_ref;
    }
    END_TIMER("outer_prod_map");
    EXPECT_DOUBLE_EQ( result_oprod_transp(1), result_oprod_map(1) );
    END_TIMER("OPROD");
    // 0.0311495 vs. 0.0308666

    this->profiler_output("arena_allocation");
}
