#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

TEST(Eigen, outer_product_scalar) {
    typedef Eigen::Array<double,Eigen::Dynamic,1>  ArrayDbl;
    typedef Eigen::Vector<ArrayDbl,Eigen::Dynamic> TableDbl;

    TableDbl data_table;
    data_table.resize(10);
    // lengths of data_table rows
    std::vector<unsigned int> row_len = {2, 4, 8};
    // declare ref values
    ArrayDbl ref_weights(4);
    ref_weights << 0.5, 0.6, 0.7, 0.8;

    // resize and fill data
    for (uint i=0; i<row_len.size(); ++i)
        data_table(i).resize(row_len[i]);
    data_table(0) << 1.1, 1.2;                  // element determinants
    data_table(1) << ref_weights;               // weights are copied from ref data
    data_table(2) << 0, 0, 0, 0, 0, 0, 0, 0;    // JxW (reset values)

    // Simulation of JxW evaluation: weights * determinant
    Eigen::Map<Eigen::Vector<double, 4>> weight_vals(data_table(1).data(), 4);   // 4 = n quadrature points
    Eigen::Map<Eigen::Vector<double, 2>> det_value(data_table(0).data(), 2);     // 2 = n elements on patch
    Eigen::Map<Eigen::Matrix<double, 4, 2>> result_val(data_table(2).data(), 4, 2);
    result_val = weight_vals * det_value.transpose();

    // Output - product
    std::cout << "Weights:\n" << weight_vals << std::endl << std::endl;
    std::cout << "Det:\n" << det_value << std::endl << std::endl;
    std::cout << "JxW:\n" << result_val << std::endl << std::endl;

    // Output - data table
    std::cout << "Output:" << std::endl;
    for (uint i=0; i<row_len.size(); ++i) {
        for (uint j=0; j<row_len[i]; ++j)
            std::cout << data_table(i)(j) << " ";
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
}


TEST(Eigen, outer_product_vector) {
    typedef Eigen::Array<double,Eigen::Dynamic,1>  ArrayDbl;
    typedef Eigen::Vector<ArrayDbl,Eigen::Dynamic> TableDbl;

    TableDbl data_table;
    data_table.resize(15);
    // lengths of data_table rows
    std::vector<unsigned int> row_len = {2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 12, 12, 24, 24, 24};
    for (uint i=0; i<row_len.size(); ++i) {
        data_table(i).resize(row_len[i]);
        // fill 'shorter' rows, 'longer' rows will be computed
        for (uint j=0; j<row_len[i]; ++j)
            if (i<12) data_table(i)(j) = i + 1.1 + j*0.1;
            else data_table(i)(j) = 0.0;
    }

    // Simulation of JxW evaluation: jac_value.transpose() * ref_value
    Eigen::Map<Eigen::Matrix<ArrayDbl, Eigen::Dynamic, Eigen::Dynamic>> jac_value(data_table.data() + 0, 3, 3);
    Eigen::Map<Eigen::Matrix<ArrayDbl, Eigen::Dynamic, Eigen::Dynamic>> ref_value(data_table.data() + 9, 3, 1);
    std::cout << "ref value 0:\n" << ref_value(0) << std::endl;
    std::cout << "ref value 1:\n" << ref_value(1) << std::endl;
    std::cout << "ref value 2:\n" << ref_value(2) << std::endl;
//    Eigen::Vector<double, 4> weight_vals;
//    weight_vals << 0.5, 0.6, 0.7, 0.8;
//    Eigen::Map<Eigen::Vector<double, 2>> det_value(data_table(0).data(), 2);
//    std::cout << "Ref vals:\n" << ref_vals << std::endl << std::endl;
//    std::cout << "Det:\n" << det_value << std::endl << std::endl;
//    Eigen::Map<Eigen::Matrix<double, 4, 2>> result_val(data_table(7).data(), 4, 2);
//    result_val = weight_vals * det_value.transpose();
//    std::cout << "Result:\n" << result_val << std::endl << std::endl;

    // Output
    std::cout << "Output:" << std::endl;
    for (uint i=0; i<row_len.size(); ++i) {
        for (uint j=0; j<row_len[i]; ++j)
            std::cout << data_table(i)(j) << " ";
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
}


//using Eigen::MatrixXd;
//using Eigen::VectorXd;
//
//TEST(Eigen, multi_simple) {
//    MatrixXd m;
//    m.resize(3,3);
//    m << 1, 2, 3, 2, 4, 5, 3, 5, 6;
//    std::cout << "m =" << std::endl << m << std::endl;
//    VectorXd v(3);
//    v << 1, 2, 3;
//    std::cout << "m * v =" << std::endl << m * v << std::endl;
//    std::cout << "det   =" << std::endl << m.determinant() << std::endl;
//    std::cout << "inv m =" << std::endl << m.inverse() << std::endl;
//}
//
//
//template <unsigned int size>
//class VectorCol {
//public:
//    typename Eigen::Matrix<double,size,1> data;
//
//    VectorCol() {}
//
//    VectorCol(int n) {}
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
//    inline VectorCol<size> operator-(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//        res.data = this->data - other.data;
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
//        for (unsigned int i=0; i<size; ++i)
//            res.data[i] = this->data[i] * other.data[i];
//        return res;
//    }
//
//    inline VectorCol<size> operator/(const VectorCol<size> &other) const {
//    	VectorCol<size> res;
//        for (unsigned int i=0; i<size; ++i)
//            res.data[i] = this->data[i] / other.data[i];
//        return res;
//    }
//
//};
//
//class Vector16 {
//public:
//    typename Eigen::Matrix<double,16,1> data;
//
//    Vector16() {}
//
//    //Vector16(int n) {}
//
//    inline double & operator[](std::size_t item) {
//        return data[item];
//    }
//
//    inline double & operator()(std::size_t item) {
//        return data[item];
//    }
//
////    inline Vector16 operator+(const Vector16 &other) const {
////    	Vector16 res;
////        res.data = this->data + other.data;
////        return res;
////    }
//
//    inline Vector16 operator-(const Vector16 &other) const {
//    	Vector16 res;
//        res.data = this->data - other.data;
//        return res;
//    }
//
//    inline Vector16 operator*(const double &coef) const {
//    	Vector16 res;
//        res.data = this->data * coef;
//        return res;
//    }
//
//    inline Vector16 operator/(const double &coef) const {
//    	Vector16 res;
//        res.data = this->data / coef;
//        return res;
//    }
//
//    inline Vector16 operator*(const Vector16 &other) const {
//    	Vector16 res;
//        for (unsigned int i=0; i<16; ++i)
//            res.data[i] = this->data[i] * other.data[i];
//        return res;
//    }
//
//    inline Vector16 operator/(const Vector16 &other) const {
//    	Vector16 res;
//        for (unsigned int i=0; i<16; ++i)
//            res.data[i] = this->data[i] / other.data[i];
//        return res;
//    }
//
//};
//
//
//Vector16 operator+(const Vector16 &A, const Vector16 &B) {
//	Vector16 res;
//    res.data = A.data + B.data;
//    return res;
//}
//
//
//
//class Vector16Plus : public Vector16 {
//public:
//	Vector16Plus() : Vector16() {}
//
//    //Vector16Plus(int n) : Vector16(n) {}
//};
//
//
//
//template class VectorCol<8>;
//template class VectorCol<200>;
//
//
//typedef Eigen::Vector<VectorCol<8>,3> Vector3x8;
//typedef Eigen::Matrix<VectorCol<8>,3,3> Matrix3x8;
//typedef Eigen::Vector<VectorCol<200>,3> Vector3x200;
//typedef Eigen::Matrix<VectorCol<200>,3,3> Matrix3x200;
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
//    template<> struct NumTraits<Vector16> : GenericNumTraits<Vector16>
//    {
//        typedef Vector16 Real;
//        typedef Vector16 NonInteger;
//        typedef Vector16 Nested;
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
//    template<> struct NumTraits<Vector16Plus> : GenericNumTraits<Vector16Plus>
//    {
//        typedef Vector16Plus Real;
//        typedef Vector16Plus NonInteger;
//        typedef Vector16Plus Nested;
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
//
//
//TEST(Eigen, own_types) {
//	// mat multi
//    Eigen::Matrix<Vector16,2,3> a;
//    a(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9;
//    a(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7;
//    a(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5;
//    a(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2;
//    a(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4;
//    a(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3;
//
//    Eigen::Matrix<Vector16,3,4> b;
//    b(0,0).data << 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7;
//    b(0,1).data << 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9;
//    b(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5;
//    b(0,3).data << 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3;
//    b(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2;
//    b(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4;
//    b(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3, 5, 1, 3, 3;
//    b(1,3).data << 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7;
//    b(2,0).data << 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5;
//    b(2,1).data << 4, 5, 7, 5, 9, 5, 1, 5, 4, 5, 5, 5, 7, 5, 1, 5;
//    b(2,2).data << 3, 2, 5, 2, 2, 3, 5, 2, 3, 2, 5, 2, 2, 3, 5, 2;
//    b(2,3).data << 1, 2, 3, 5, 4, 7, 6, 8, 1, 2, 3, 5, 4, 7, 6, 8;
//
//    Eigen::Matrix<Vector16,2,4> c = a * b;
//    std::cout << "c(0,0) =" << std::endl << c(0,0).data << std::endl;
//
//    //assignment operator
//    Eigen::Matrix<Vector16,2,3> x;
//    x(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9, 1, 2, 1, 9;
//    x(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7, 2, 4, 4, 7;
//    x(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5, 3, 5, 2, 5;
//    x(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5, 2;
//    x(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4, 4, 2, 1, 4;
//    x(2,2).data << 1, 2, 3, 5, 4, 7, 6, 8, 1, 2, 3, 5, 4, 7, 6, 8;
//
//    Eigen::Matrix<Vector16,2,3> y = x + a;
//    std::cout << "y(0,0) =" << std::endl << y(0,0).data << std::endl;
//
//}


//TEST(Eigen, multi_vectorized_8) {
//    Vector3x8 vec;
//    vec(0).data << 1, 1, 1, 1, 2, 2, 2, 2;
//    vec(1).data << 2, 2, 2, 2, 3, 3, 3, 3;
//    vec(2).data << 3, 3, 3, 3, 4, 4, 4, 4;
//
//    Matrix3x8 mat;
//    mat(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9;
//    mat(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7;
//    mat(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5;
//    mat(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2;
//    mat(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4;
//    mat(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3;
//    mat(2,0).data << 3, 0, 2, 1, 3, 0, 2, 1;
//    mat(2,1).data << 5, 6, 3, 3, 5, 6, 3, 3;
//    mat(2,2).data << 6, 4, 4, 8, 6, 4, 4, 8;
//
//    Vector3x8 res = mat * vec;
//    for (uint i=0; i<3; ++i) {
//        std::cout << "res[" << i << "] =" << std::endl << res[i].data << std::endl;
//    }
//
//    VectorCol<8> det = mat.determinant();
//    std::cout << "det =" << std::endl << det.data << std::endl;
//
//    Matrix3x8 inv = mat.inverse();
//    std::cout << "inv =" << std::endl << inv(0,0).data << std::endl;
//}
//
//class EigenTest : public testing::Test {
//public:
//	EigenTest()
//    {
//		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/system";
//        FilePath::set_io_dirs(".",root_dir,"",".");
//        Profiler::instance();
//        Profiler::set_memory_monitoring(false, false);
//    }
//
//    ~EigenTest()
//    {
//        Profiler::uninitialize();
//    }
//
//	/// Perform profiler output.
//    void profiler_output(std::string file_name) {
//		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
//		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
//	}
//};
//
//
//TEST_F(EigenTest, multi_vectorized_200) {
//	static const uint N_RUNS = 4e7 / 200;
//
//    Vector3x8 vec_data;
//    vec_data(0).data << 1, 1, 1, 1, 2, 2, 2, 2;
//    vec_data(1).data << 2, 2, 2, 2, 3, 3, 3, 3;
//    vec_data(2).data << 3, 3, 3, 3, 4, 4, 4, 4;
//
//    Matrix3x8 mat_data;
//    mat_data(0,0).data << 1, 2, 1, 9, 1, 2, 1, 9;
//    mat_data(0,1).data << 2, 4, 4, 7, 2, 4, 4, 7;
//    mat_data(0,2).data << 3, 5, 2, 5, 3, 5, 2, 5;
//    mat_data(1,0).data << 2, 3, 5, 2, 2, 3, 5, 2;
//    mat_data(1,1).data << 4, 2, 1, 4, 4, 2, 1, 4;
//    mat_data(1,2).data << 5, 1, 3, 3, 5, 1, 3, 3;
//    mat_data(2,0).data << 3, 0, 2, 1, 3, 0, 2, 1;
//    mat_data(2,1).data << 5, 6, 3, 3, 5, 6, 3, 3;
//    mat_data(2,2).data << 6, 4, 4, 8, 6, 4, 4, 8;
//
//    Vector3x200 vec;
//    Matrix3x200 mat;
//    for (uint i=0; i<200; ++i) {
//        uint data_i = i%8;
//    	mat(0,0).data(i) = mat_data(0,0).data(data_i);
//    	mat(0,1).data(i) = mat_data(0,1).data(data_i);
//    	mat(0,2).data(i) = mat_data(0,2).data(data_i);
//    	mat(1,0).data(i) = mat_data(1,0).data(data_i);
//    	mat(1,1).data(i) = mat_data(1,1).data(data_i);
//    	mat(1,2).data(i) = mat_data(1,2).data(data_i);
//    	mat(2,0).data(i) = mat_data(2,0).data(data_i);
//    	mat(2,1).data(i) = mat_data(2,1).data(data_i);
//    	mat(2,2).data(i) = mat_data(2,2).data(data_i);
//        vec(0).data(i)  = vec_data(0).data(data_i);
//        vec(1).data(i)  = vec_data(1).data(data_i);
//        vec(2).data(i)  = vec_data(2).data(data_i);
//    }
//
//    VectorCol<200> result_multi;
//    VectorCol<200> result_det;
//    VectorCol<200> result_inv;
//
//    START_TIMER("multi_mat_vec");
//    for (unsigned int i=0; i<N_RUNS; ++i) {
//        Vector3x200 multi = mat * vec;
//        result_multi.data = result_multi.data + multi[0].data;
//    }
//    END_TIMER("multi_mat_vec");
//
//    START_TIMER("determinant");
//    for (unsigned int i=0; i<N_RUNS; ++i) {
//        VectorCol<200> det = mat.determinant();
//        result_det.data = result_det.data + det.data;
//    }
//    END_TIMER("determinant");
//
//    START_TIMER("inverse");
//    for (unsigned int i=0; i<N_RUNS; ++i) {
//        Matrix3x200 inv = mat.inverse();
//        result_multi.data = result_multi.data + inv(0,0).data + inv(1,1).data + inv(2,2).data;
//    }
//    END_TIMER("inverse");
//
//    this->profiler_output("eigen");
//}
