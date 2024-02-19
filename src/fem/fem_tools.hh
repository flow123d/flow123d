/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    mapping.hh
 * @brief   Class Mapping calculates data related to the mapping
 *          of the reference cell to the actual cell, such as Jacobian
 *          and normal vectors.
 * @author  Jan Stebel
 */

#ifndef MAPPING_HH_
#define MAPPING_HH_

#include <armadillo>
#include <vector>
#include <limits>
#include "system/fmt/posix.h"           // for FMT_UNUSED




//#include <Eigen/Core>
//#include <Eigen/Dense>
//
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


namespace fe_tools {



/**
 * @brief Calculates determinant of a rectangular matrix.
 */
template<class T>
double determinant(const T &M);



inline arma::mat::fixed<1,1> normal_matrix(const arma::mat::fixed<1,2> &A) {
    arma::mat::fixed<1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1);
    return res;
}

inline arma::mat::fixed<1,1> normal_matrix(const arma::mat::fixed<2,1> &A) {
    arma::mat::fixed<1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0);
    return res;
}

inline arma::mat::fixed<1,1> normal_matrix(const arma::mat::fixed<1,3> &A) {
    arma::mat::fixed<1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    return res;
}

inline arma::mat::fixed<1,1> normal_matrix(const arma::mat::fixed<3,1> &A) {
    arma::mat::fixed<1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    return res;
}

inline arma::mat::fixed<2,2> normal_matrix(const arma::mat::fixed<2,3> &A) {
    arma::mat::fixed<2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    res(0,1) = A(0,0)*A(1,0)+A(0,1)*A(1,1)+A(0,2)*A(1,2);
    res(1,0) = A(1,0)*A(0,0)+A(1,1)*A(0,1)+A(1,2)*A(0,2);
    res(1,1) = A(1,0)*A(1,0)+A(1,1)*A(1,1)+A(1,2)*A(1,2);
    return res;
}

inline arma::mat::fixed<2,2> normal_matrix(const arma::mat::fixed<3,2> &A) {
    arma::mat::fixed<2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    res(0,1) = A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    res(1,0) = A(0,1)*A(0,0)+A(1,1)*A(1,0)+A(2,1)*A(2,0);
    res(1,1) = A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    return res;
}



template<> inline double determinant(const arma::mat::fixed<1,1> &M)
{
    return M(0,0);
}

template<> inline double determinant(const arma::mat::fixed<2,2> &M)
{
    return M(0,0)*M(1,1) - M(1,0)*M(0,1);
}

template<> inline double determinant(const arma::mat::fixed<3,3> &M)
{
    return ( M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) )
         - ( M(2,0)*M(1,1)*M(0,2) + M(2,1)*M(1,2)*M(0,0) + M(2,2)*M(1,0)*M(0,1) );
}

template<> inline double determinant(FMT_UNUSED const arma::mat::fixed<0,3> &M)
{
    return 0;
}

template<> inline double determinant(FMT_UNUSED const arma::mat::fixed<3,0> &M)
{
    return 0;
}

template<> inline double determinant(const arma::mat::fixed<1,2> &M)
{
    return sqrt( determinant(normal_matrix(M)) );
}

template<> inline double determinant(const arma::mat::fixed<2,1> &M)
{
    return sqrt( determinant(normal_matrix(M)) );
}

template<> inline double determinant(const arma::mat::fixed<1,3> &M)
{
    return sqrt( determinant(normal_matrix(M)) );
}

template<> inline double determinant(const arma::mat::fixed<3,1> &M)
{
    return sqrt( determinant(normal_matrix(M)) );
}

template<> inline double determinant(const arma::mat::fixed<2,3> &M)
{
    return sqrt( determinant(normal_matrix(M)) );
}

template<> inline double determinant(const arma::mat::fixed<3,2> &M)
{
	return sqrt( determinant(normal_matrix(M)) );
}


/**
 * @brief Calculates inverse of rectangular matrix or pseudoinverse of non-rectangular matrix.
 */
template<arma::uword m, arma::uword n>
arma::mat::fixed<n,m> inverse(const arma::mat::fixed<m,n> &A) {
    if (m<n) return A.t() * inverse(normal_matrix(A));
    else return inverse(normal_matrix(A)) * A.t();
}


template<> inline arma::mat::fixed<1,1> inverse(const arma::mat::fixed<1,1> &A)
{
	arma::mat::fixed<1,1> B;
	B(0,0) = 1/A(0,0);
	return B;
}

template<> inline arma::mat::fixed<2,2> inverse(const arma::mat::fixed<2,2> &A)
{
	arma::mat::fixed<2,2> B;
	double det = determinant(A);

	B(0,0) = A(1,1) / det;
	B(0,1) = -A(0,1) / det;
	B(1,0) = -A(1,0) / det;
	B(1,1) = A(0,0) / det;
	return B;
}

template<> inline arma::mat::fixed<3,3> inverse(const arma::mat::fixed<3,3> &A)
{
	arma::mat::fixed<3,3> B;

	B(0,0) = A(1,1)*A(2,2) - A(2,1)*A(1,2);
	B(1,0) = -A(1,0)*A(2,2) + A(2,0)*A(1,2);
	B(2,0) = A(1,0)*A(2,1) - A(2,0)*A(1,1);

	double det = A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0);
	B(0,0)/=det; B(1,0)/=det; B(2,0)/=det;

	B(0,1) = (-A(0,1)*A(2,2) + A(2,1)*A(0,2)) / det;
	B(1,1) = (A(0,0)*A(2,2) - A(2,0)*A(0,2)) / det;
	B(2,1) = (-A(0,0)*A(2,1) + A(2,0)*A(0,1)) / det;

	B(0,2) = (A(0,1)*A(1,2) - A(1,1)*A(0,2)) / det;
	B(1,2) = (-A(0,0)*A(1,2) + A(1,0)*A(0,2)) / det;
	B(2,2) = (A(0,0)*A(1,1) - A(1,0)*A(0,1)) / det;

	return B;
}


} // closing namespace fe_tools



//namespace eigen_tools {
//
//// use only one of the following typedef
//typedef VectorCol<200> Vec200;
////typedef Eigen::Array<double,200,1> Vec200;
//
//
//
///**
// * @brief Calculates determinant of a rectangular matrix.
// */
//template<class T>
//Vec200 determinant(const T &M);
//
//
//
//inline Eigen::Matrix<Vec200,1,1> normal_matrix(const Eigen::Matrix<Vec200,1,2> &A) {
//	Eigen::Matrix<Vec200,1,1> res;
//    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1);
//    return res;
//}
//
//inline Eigen::Matrix<Vec200,1,1> normal_matrix(const Eigen::Matrix<Vec200,2,1> &A) {
//	Eigen::Matrix<Vec200,1,1> res;
//    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0);
//    return res;
//}
//
//inline Eigen::Matrix<Vec200,1,1> normal_matrix(const Eigen::Matrix<Vec200,1,3> &A) {
//	Eigen::Matrix<Vec200,1,1> res;
//    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
//    return res;
//}
//
//inline Eigen::Matrix<Vec200,1,1> normal_matrix(const Eigen::Matrix<Vec200,3,1> &A) {
//	Eigen::Matrix<Vec200,1,1> res;
//    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
//    return res;
//}
//
//inline Eigen::Matrix<Vec200,2,2> normal_matrix(const Eigen::Matrix<Vec200,2,3> &A) {
//    Eigen::Matrix<Vec200,2,2> res;
//    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
//    res(0,1) = A(0,0)*A(1,0)+A(0,1)*A(1,1)+A(0,2)*A(1,2);
//    res(1,0) = A(1,0)*A(0,0)+A(1,1)*A(0,1)+A(1,2)*A(0,2);
//    res(1,1) = A(1,0)*A(1,0)+A(1,1)*A(1,1)+A(1,2)*A(1,2);
//    return res;
//}
//
//inline Eigen::Matrix<Vec200,2,2> normal_matrix(const Eigen::Matrix<Vec200,3,2> &A) {
//	Eigen::Matrix<Vec200,2,2> res;
//    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
//    res(0,1) = A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
//    res(1,0) = A(0,1)*A(0,0)+A(1,1)*A(1,0)+A(2,1)*A(2,0);
//    res(1,1) = A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
//    return res;
//}
//
//
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,1,1> &M)
//{
//    return M(0,0);
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,2,2> &M)
//{
//    return M(0,0)*M(1,1) - M(1,0)*M(0,1);
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,3,3> &M)
//{
//    return ( M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) )
//         - ( M(2,0)*M(1,1)*M(0,2) + M(2,1)*M(1,2)*M(0,0) + M(2,2)*M(1,0)*M(0,1) );
//}
//
//template<> inline Vec200 determinant(FMT_UNUSED const Eigen::Matrix<Vec200,0,3> &M)
//{
//    return Vec200();
//}
//
//template<> inline Vec200 determinant(FMT_UNUSED const Eigen::Matrix<Vec200,3,0> &M)
//{
//    return Vec200();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,1,2> &M)
//{
//    return determinant(normal_matrix(M)).sqrt();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,2,1> &M)
//{
//    return determinant(normal_matrix(M)).sqrt();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,1,3> &M)
//{
//    return determinant(normal_matrix(M)).sqrt();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,3,1> &M)
//{
//    return determinant(normal_matrix(M)).sqrt();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,2,3> &M)
//{
//    return determinant(normal_matrix(M)).sqrt();
//}
//
//template<> inline Vec200 determinant(const Eigen::Matrix<Vec200,3,2> &M)
//{
//	return determinant(normal_matrix(M)).sqrt();
//}
//
//
///**
// * @brief Calculates inverse of rectangular matrix or pseudoinverse of non-rectangular matrix.
// */
//template<int m, int n>
//Eigen::Matrix<Vec200,n,m> inverse(const Eigen::Matrix<Vec200,m,n> &A) {
//    // only for cases m > n
//    return inverse(normal_matrix(A)) * A.transpose();
//}
//
//
//template<> inline Eigen::Matrix<Vec200,1,1> inverse<1,1>(const Eigen::Matrix<Vec200,1,1> &A)
//{
//	Eigen::Matrix<Vec200,1,1> B;
//    B(0,0) = A(0,0).inverse(); // 1/A(0,0)
//    return B;
//}
//
//template<> inline Eigen::Matrix<Vec200,2,2> inverse<2,2>(const Eigen::Matrix<Vec200,2,2> &A)
//{
//	Eigen::Matrix<Vec200,2,2> B;
//	Vec200 det = determinant(A);
//
//    B(0,0) = A(1,1) / det;
//    B(0,1) = -A(0,1) / det;
//    B(1,0) = -A(1,0) / det;
//    B(1,1) = A(0,0) / det;
//    return B;
//}
//
//template<> inline Eigen::Matrix<Vec200,3,3> inverse<3,3>(const Eigen::Matrix<Vec200,3,3> &A)
//{
//    Eigen::Matrix<Vec200,3,3> B;
//
//    B(0,0) = A(1,1)*A(2,2) - A(2,1)*A(1,2);
//    B(0,1) = A(2,0)*A(1,2) - A(1,0)*A(2,2);
//    B(0,2) = A(1,0)*A(2,1) - A(2,0)*A(1,1);
//
//    Vec200 det = A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
//    B(0,0) = B(0,0) / det;
//    B(0,1) = B(0,1) / det;
//    B(0,2) = B(0,2) / det;
//
//    B(1,0) = (A(2,1)*A(0,2) - A(0,1)*A(2,2)) / det;
//    B(1,1) = (A(0,0)*A(2,2) - A(2,0)*A(0,2)) / det;
//    B(1,2) = (A(2,0)*A(0,1) - A(0,0)*A(2,1)) / det;
//
//    B(2,0) = (A(0,1)*A(1,2) - A(1,1)*A(0,2)) / det;
//    B(2,1) = (A(1,0)*A(0,2) - A(0,0)*A(1,2)) / det;
//    B(2,2) = (A(0,0)*A(1,1) - A(1,0)*A(0,1)) / det;
//
//    return B;
//}
//
//template<> inline Eigen::Matrix<Vec200,2,1> inverse<1,2>(const Eigen::Matrix<Vec200,1,2> &A)
//{
//    return A.transpose() * inverse(normal_matrix(A));
//}
//
//template<> inline Eigen::Matrix<Vec200,3,1> inverse<1,3>(const Eigen::Matrix<Vec200,1,3> &A)
//{
//    return A.transpose() * inverse(normal_matrix(A));
//}
//
//template<> inline Eigen::Matrix<Vec200,3,2> inverse<2,3>(const Eigen::Matrix<Vec200,2,3> &A)
//{
//    return A.transpose() * inverse(normal_matrix(A));
//}
//
//
//} // closing namespace eigen_tools



#endif /* MAPPING_HH_ */
