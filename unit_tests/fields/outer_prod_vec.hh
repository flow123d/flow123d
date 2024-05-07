#ifndef OUTER_PRODUCT_HH_
#define OUTER_PRODUCT_HH_

#include <armadillo>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "system/asserts.hh"


class OuterProdVec : public Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
{
public:
    using Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Matrix; // Inherit constructors

    // Override the operator+
    OuterProdVec operator+(const OuterProdVec& other) const
    {
        ASSERT_EQ(this->rows(), other.rows());
        ASSERT_EQ(this->cols(), other.cols());
    	return OuterProdVec(static_cast<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&>(*this) + other);
    }

    // Override the operator*
    OuterProdVec operator*(const OuterProdVec& other) const
    {
        ASSERT_EQ(this->cols(), 1);
        ASSERT_EQ(other.cols(), 1);
        return OuterProdVec(static_cast<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&>(*this) * other.transpose());
    }
};


int test_add_operator_simple()
{
    OuterProdVec arr1(3, 1);
    OuterProdVec arr2(3, 1);

    // Initialize arrays with some values
    arr1 << 1, 2, 3;
    arr2 << 4, 5, 6;

    // Add arrays using overridden operator+
    OuterProdVec result = arr1 + arr2;

    // Output the result
    std::cout << "Result:\n" << result << std::endl;

    return 0;
}

int test_multi_operator_simple()
{
    OuterProdVec arr1(3, 1);
    OuterProdVec arr2(3, 1);

    // Initialize arrays with some values
    arr1 << 1, 2, 3;
    arr2 << 4, 5, 6;

    // Multiple arrays using overridden operator*
    OuterProdVec result = arr1 * arr2;

    // Output the result
    std::cout << "Result:\n" << result << std::endl;

    return 0;
}

int test_add_operator()
{
	Eigen::Vector<OuterProdVec, 3> arr1;
	Eigen::Vector<OuterProdVec, 3> arr2;

    // Initialize arrays with some values
	for (uint i=0; i<3; ++i) {
        arr1(i).resize(3, 1);
        arr1(i) << (i+1), (i+2), (i+3);
        arr2(i).resize(3, 1);
        arr2(i) << (2*i+1), (2*i+2), (2*i+3);
	}

    // Add arrays using overridden operator+
    Eigen::Vector<OuterProdVec, 3> result = arr1 + arr2;

    // Output the result
    std::cout << "Result 0:\n" << result(0) << std::endl;
    std::cout << "Result 1:\n" << result(1) << std::endl;
    std::cout << "Result 2:\n" << result(2) << std::endl;

    return 0;
}

int test_multi_operator()
{
    Eigen::Array<OuterProdVec, 3, 1> arr1;
    Eigen::Array<OuterProdVec, 3, 1> arr2;

    // Initialize arrays with some values
    for (uint i=0; i<3; ++i) {
        arr1(i).resize(3, 1);
        arr1(i) << (i+1), (i+2), (i+3);
        arr2(i).resize(3, 1);
        arr2(i) << (2*i+1), (2*i+2), (2*i+3);
    }

    // Multiple arrays using overridden operator*
    //auto result = arr1 * arr2; // ERROR, during compilation, see log

    // compilable and functional code, fix of line 107
    Eigen::Array<OuterProdVec, 3, 3> result;
    for (uint i=0; i<result.rows(); ++i)
        for (uint j=0; j<result.cols(); ++j)
            result(i,j) = arr1(i) * arr2(j);

    // Output the result
    std::cout << "Array1 0:\n" << arr1(0) << std::endl;
    std::cout << "Array1 1:\n" << arr1(1) << std::endl;
    std::cout << "Array2 0:\n" << arr2(0) << std::endl;
    std::cout << "Array2 1:\n" << arr2(1) << std::endl;
    std::cout << "Result 0:\n" << result(0,0) << std::endl;
    std::cout << "Result 1:\n" << result(1,1) << std::endl;

    return 0;
}

#endif /* OUTER_PRODUCT_HH_ */
