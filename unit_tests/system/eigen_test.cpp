#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "system/file_path.hh"
#include "system/sys_profiler.hh"


using Eigen::MatrixXd;
using Eigen::VectorXd;

//TEST(Eigen, multi_simple) {
//    MatrixXd m(3,3);
//    m << 1, 2, 3, 2, 4, 5, 3, 5, 6;
//    std::cout << "m =" << std::endl << m << std::endl;
//    VectorXd v(3);
//    v << 1, 2, 3;
//    std::cout << "m * v =" << std::endl << m * v << std::endl;
//    std::cout << "det   =" << std::endl << m.determinant() << std::endl;
//    std::cout << "inv m =" << std::endl << m.inverse() << std::endl;
//}


template <unsigned int size>
class VectorCol {
public:
    typename Eigen::Matrix<double,size,1> data;

    VectorCol() {}

    VectorCol(int n) {}

    inline double & operator[](std::size_t item) {
        return data[item];
    }

    inline double & operator()(std::size_t item) {
        return data[item];
    }

    inline VectorCol<size> operator+(const VectorCol<size> &other) const {
    	VectorCol<size> res;
        res.data = this->data + other.data;
        return res;
    }

    inline VectorCol<size> operator-(const VectorCol<size> &other) const {
    	VectorCol<size> res;
        res.data = this->data - other.data;
        return res;
    }

    inline VectorCol<size> operator*(const double &coef) const {
    	VectorCol<size> res;
        res.data = this->data * coef;
        return res;
    }

    inline VectorCol<size> operator/(const double &coef) const {
    	VectorCol<size> res;
        res.data = this->data / coef;
        return res;
    }

    inline VectorCol<size> operator*(const VectorCol<size> &other) const {
    	VectorCol<size> res;
        for (unsigned int i=0; i<size; ++i)
            res.data[i] = this->data[i] * other.data[i];
        return res;
    }

    inline VectorCol<size> operator/(const VectorCol<size> &other) const {
    	VectorCol<size> res;
        for (unsigned int i=0; i<size; ++i)
            res.data[i] = this->data[i] / other.data[i];
        return res;
    }

};

template class VectorCol<8>;
template class VectorCol<200>;


typedef Eigen::Matrix<VectorCol<8>,3,1> Vector3x8;
typedef Eigen::Matrix<VectorCol<8>,3,3> Matrix3x8;
typedef Eigen::Matrix<VectorCol<200>,3,1> Vector3x200;
typedef Eigen::Matrix<VectorCol<200>,3,3> Matrix3x200;

namespace Eigen {
    template<> struct NumTraits<VectorCol<8>> : GenericNumTraits<VectorCol<8>>
    {
        typedef VectorCol<8> Real;
        typedef VectorCol<8> NonInteger;
        typedef VectorCol<8> Nested;

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 150,
            MulCost = 100
        };
    };
    template<> struct NumTraits<VectorCol<200>> : GenericNumTraits<VectorCol<200>>
    {
        typedef VectorCol<200> Real;
        typedef VectorCol<200> NonInteger;
        typedef VectorCol<200> Nested;

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = 6,
            AddCost = 150,
            MulCost = 100
        };
    };
}



TEST(Eigen, multi_vectorized_8) {
    Vector3x8 vec;
    vec(0).data << 1, 1, 1, 1, 2, 2, 2, 2;
    vec(1).data << 2, 2, 2, 2, 3, 3, 3, 3;
    vec(2).data << 3, 3, 3, 3, 4, 4, 4, 4;

    Matrix3x8 mat;
    mat(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9;
    mat(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7;
    mat(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5;
    mat(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2;
    mat(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4;
    mat(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3;
    mat(2,0).data << 3, 0, 2, 1, 3, 0, 2, 1;
    mat(2,1).data << 5, 6, 3, 3, 5, 6, 3, 3;
    mat(2,2).data << 6, 4, 4, 8, 6, 4, 4, 8;

    Vector3x8 res = mat * vec;
    for (uint i=0; i<3; ++i) {
        std::cout << "res[" << i << "] =" << std::endl << res[i].data << std::endl;
    }

    VectorCol<8> det = mat.determinant();
    std::cout << "det =" << std::endl << det.data << std::endl;

    Matrix3x8 inv = mat.inverse();
    std::cout << "inv =" << std::endl << inv(0,0).data << std::endl;
}

class EigenTest : public testing::Test {
public:
	EigenTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/system";
        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::set_memory_monitoring(false, false);
    }

    ~EigenTest()
    {
        Profiler::uninitialize();
    }

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};


TEST_F(EigenTest, multi_vectorized_200) {
	static const uint N_RUNS = 4e7 / 200;

    Vector3x8 vec_data;
    vec_data(0).data << 1, 1, 1, 1, 2, 2, 2, 2;
    vec_data(1).data << 2, 2, 2, 2, 3, 3, 3, 3;
    vec_data(2).data << 3, 3, 3, 3, 4, 4, 4, 4;

    Matrix3x8 mat_data;
    mat_data(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9;
    mat_data(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7;
    mat_data(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5;
    mat_data(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2;
    mat_data(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4;
    mat_data(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3;
    mat_data(2,0).data << 3, 0, 2, 1, 3, 0, 2, 1;
    mat_data(2,1).data << 5, 6, 3, 3, 5, 6, 3, 3;
    mat_data(2,2).data << 6, 4, 4, 8, 6, 4, 4, 8;

    Vector3x200 vec;
    Matrix3x200 mat;
    for (uint i=0; i<200; ++i) {
        uint data_i = i%8;
    	mat(0,0).data(i) = mat_data(0,0).data(data_i);
    	mat(0,1).data(i) = mat_data(0,1).data(data_i);
    	mat(0,2).data(i) = mat_data(0,2).data(data_i);
    	mat(1,0).data(i) = mat_data(1,0).data(data_i);
    	mat(1,1).data(i) = mat_data(1,1).data(data_i);
    	mat(1,2).data(i) = mat_data(1,2).data(data_i);
    	mat(2,0).data(i) = mat_data(2,0).data(data_i);
    	mat(2,1).data(i) = mat_data(2,1).data(data_i);
    	mat(2,2).data(i) = mat_data(2,2).data(data_i);
        vec(0).data(i)  = vec_data(0).data(data_i);
        vec(1).data(i)  = vec_data(1).data(data_i);
        vec(2).data(i)  = vec_data(2).data(data_i);
    }

    VectorCol<200> result_multi;
    VectorCol<200> result_det;
    VectorCol<200> result_inv;

    START_TIMER("multi_mat_vec");
    for (unsigned int i=0; i<N_RUNS; ++i) {
        Vector3x200 multi = mat * vec;
        result_multi.data = result_multi.data + multi[0].data;
    }
    END_TIMER("multi_mat_vec");

    START_TIMER("determinant");
    for (unsigned int i=0; i<N_RUNS; ++i) {
        VectorCol<200> det = mat.determinant();
        result_det.data = result_det.data + det.data;
    }
    END_TIMER("determinant");

    START_TIMER("inverse");
    for (unsigned int i=0; i<N_RUNS; ++i) {
        Matrix3x200 inv = mat.inverse();
        result_multi.data = result_multi.data + inv(0,0).data + inv(1,1).data + inv(2,2).data;
    }
    END_TIMER("inverse");

    this->profiler_output("eigen");
}
