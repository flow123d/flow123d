#ifndef OUTER_PRODUCT_HH_
#define OUTER_PRODUCT_HH_

#include <armadillo>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "system/asserts.hh"


class OuterProdVec
{
public:
    typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data;

    OuterProdVec() {}

    OuterProdVec(int n) {
        data.resize(n, 1);
    }

    OuterProdVec(int m, int n) {
        data.resize(m, n);
    }

    inline double & operator()(std::size_t i, std::size_t j) {
        return data(i, j);
    }

    // Override the operator+
    inline OuterProdVec operator+(const OuterProdVec& other) const
    {
        ASSERT_EQ(this->data.rows(), other.data.rows());
        ASSERT_EQ(this->data.cols(), other.data.cols());
        OuterProdVec res;
        res.data = this->data + other.data;
        return res;
    }

    // Override the operator*
    inline OuterProdVec operator*(const OuterProdVec& other) const
    {
//        ASSERT_EQ(this->data.cols(), 1);
//        ASSERT_EQ(other.data.cols(), 1);
        OuterProdVec res;
        res.data = this->data * other.data;
        return res;
    }
};


int test_add_operator_simple()
{
    std::cout << "test_add_operator_simple" << std::endl;

    OuterProdVec arr1(3, 1);
    OuterProdVec arr2(3, 1);

    // Initialize arrays with some values
    arr1.data << 1, 2, 3;
    arr2.data << 4, 5, 6;

    // Add arrays using overridden operator+
    OuterProdVec result = arr1 + arr2;

    // Output the result
    std::cout << "Result:\n" << result.data << std::endl;

    return 0;
}

int test_multi_operator_simple()
{
    std::cout << "test_multi_operator_simple" << std::endl;

    OuterProdVec arr1(3, 1);
    OuterProdVec arr2(1, 3);

    // Initialize arrays with some values
    arr1.data << 1, 2, 3;
    arr2.data << 4, 5, 6;

    // Multiple arrays using overridden operator*
    OuterProdVec result = arr1 * arr2;

    // Output the result
    std::cout << "Result:\n" << result.data << std::endl;

    return 0;
}

int test_add_operator()
{
    std::cout << "test_add_operator" << std::endl;

	Eigen::Vector<OuterProdVec, 3> arr1;
	Eigen::Vector<OuterProdVec, 3> arr2;

    // Initialize arrays with some values
	for (uint i=0; i<3; ++i) {
        arr1(i).data.resize(3, 1);
        arr1(i).data << (i+1), (i+2), (i+3);
        arr2(i).data.resize(3, 1);
        arr2(i).data << (2*i+1), (2*i+2), (2*i+3);
	}

    // Add arrays using overridden operator+
    Eigen::Vector<OuterProdVec, 3> result = arr1 + arr2;

    // Output the result
    std::cout << "Result 0:\n" << result(0).data << std::endl;
    std::cout << "Result 1:\n" << result(1).data << std::endl;
    std::cout << "Result 2:\n" << result(2).data << std::endl;

    return 0;
}

int test_multi_operator()
{
    std::cout << "test_multi_operator" << std::endl;

    Eigen::Array<OuterProdVec, 3, 1> arr1;
    Eigen::Array<OuterProdVec, 3, 1> arr2;

    // Initialize arrays with some values
    for (uint i=0; i<3; ++i) {
        arr1(i).data.resize(3, 1);
        arr1(i).data << (i+1), (i+2), (i+3);
        arr2(i).data.resize(1, 3);
        arr2(i).data << (2*i+1), (2*i+2), (2*i+3);
    }

    // Multiple arrays using overridden operator*
    auto result = arr1 * arr2;

    // Output the result
    std::cout << "Array1 0:\n" << arr1(0).data << std::endl;
    std::cout << "Array1 1:\n" << arr1(1).data << std::endl;
    std::cout << "Array1 2:\n" << arr1(2).data << std::endl;
    std::cout << "Array2 0:\n" << arr2(0).data << std::endl;
    std::cout << "Array2 1:\n" << arr2(1).data << std::endl;
    std::cout << "Array2 2:\n" << arr2(2).data << std::endl;
    std::cout << "Result 0:\n" << result(0,0).data << std::endl;
    std::cout << "Result 1:\n" << result(1,0).data << std::endl;
    std::cout << "Result 2:\n" << result(2,0).data << std::endl;
    std::cout << "Result " << result.rows() << " - " << result.cols() << std::endl;

    return 0;
}

int test_multi_mat_vec_operator()
{
    std::cout << "test_multi_mat_vec_operator" << std::endl;

    Eigen::Matrix<OuterProdVec, 3, 3> arr1;
    Eigen::Matrix<OuterProdVec, 3, 1> arr2;

    // Initialize arrays with some values
    for (uint i=0; i<3; ++i) {
    	for (uint j=0; j<3; ++j) {
            arr1(i,j).data.resize(3, 3);
            arr1(i,j).data << (i+j+1), (i+j+2), (i+j+3),
                              (i+j+3), (i+j+4), (i+j+5),
		                      (i+j+5), (i+j+6), (i+j+7)/* */;
    	}
        arr2(i).data.resize(3, 1);
        arr2(i).data << (2*i+1), (2*i+2), (2*i+3);
    }

    // Multiple arrays using overridden operator*
    auto result = arr1 * arr2;

    // Output the result
    std::cout << "Array1 0,0:\n" << arr1(0,0).data << std::endl;
    std::cout << "Array1 1,0:\n" << arr1(1,0).data << std::endl;
    std::cout << "Array1 2,0:\n" << arr1(2,0).data << std::endl;
    std::cout << "Array1 0,1:\n" << arr1(0,1).data << std::endl;
    std::cout << "Array1 1,1:\n" << arr1(1,1).data << std::endl;
    std::cout << "Array1 2,1:\n" << arr1(2,1).data << std::endl;
    std::cout << "Array1 0,2:\n" << arr1(0,2).data << std::endl;
    std::cout << "Array1 1,2:\n" << arr1(1,2).data << std::endl;
    std::cout << "Array1 2,2:\n" << arr1(2,2).data << std::endl;
    std::cout << "Array2 0:\n" << arr2(0).data << std::endl;
    std::cout << "Array2 1:\n" << arr2(1).data << std::endl;
    std::cout << "Array2 2:\n" << arr2(2).data << std::endl;
    std::cout << "Result:" << result.rows() << ", " << result.cols() << std::endl;
    std::cout << "Result 0:\n" << result(0).data << std::endl;
    std::cout << "Result 1:\n" << result(1).data << std::endl;
    std::cout << "Result 2:\n" << result(2).data << std::endl;

    return 0;
}

#endif /* OUTER_PRODUCT_HH_ */
