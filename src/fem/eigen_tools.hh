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
 * @file    eigen_tools.hh
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#ifndef PATCH_DATA_TABLE_HH_
#define PATCH_DATA_TABLE_HH_

#include <Eigen/Core>
#include <Eigen/Dense>


/// Definitions of Eigen structures
typedef Eigen::Array<double,300,1>             ArrayDbl;
typedef Eigen::Array<int,300,1>                ArrayInt;
typedef Eigen::Vector<ArrayDbl,Eigen::Dynamic> TableDbl;
typedef Eigen::Vector<ArrayInt,Eigen::Dynamic> TableInt;


namespace eigen_tools {


/**
 * @brief Calculates determinant of a rectangular matrix.
 */
template<class T>
ArrayDbl determinant(const T &M);



inline Eigen::Matrix<ArrayDbl,1,1> normal_matrix(const Eigen::Matrix<ArrayDbl,1,2> &A) {
	Eigen::Matrix<ArrayDbl,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1);
    return res;
}

inline Eigen::Matrix<ArrayDbl,1,1> normal_matrix(const Eigen::Matrix<ArrayDbl,2,1> &A) {
	Eigen::Matrix<ArrayDbl,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0);
    return res;
}

inline Eigen::Matrix<ArrayDbl,1,1> normal_matrix(const Eigen::Matrix<ArrayDbl,1,3> &A) {
	Eigen::Matrix<ArrayDbl,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    return res;
}

inline Eigen::Matrix<ArrayDbl,1,1> normal_matrix(const Eigen::Matrix<ArrayDbl,3,1> &A) {
	Eigen::Matrix<ArrayDbl,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    return res;
}

inline Eigen::Matrix<ArrayDbl,2,2> normal_matrix(const Eigen::Matrix<ArrayDbl,2,3> &A) {
    Eigen::Matrix<ArrayDbl,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    res(0,1) = A(0,0)*A(1,0)+A(0,1)*A(1,1)+A(0,2)*A(1,2);
    res(1,0) = A(1,0)*A(0,0)+A(1,1)*A(0,1)+A(1,2)*A(0,2);
    res(1,1) = A(1,0)*A(1,0)+A(1,1)*A(1,1)+A(1,2)*A(1,2);
    return res;
}

inline Eigen::Matrix<ArrayDbl,2,2> normal_matrix(const Eigen::Matrix<ArrayDbl,3,2> &A) {
	Eigen::Matrix<ArrayDbl,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    res(0,1) = A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    res(1,0) = A(0,1)*A(0,0)+A(1,1)*A(1,0)+A(2,1)*A(2,0);
    res(1,1) = A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    return res;
}



template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,1,1> &M)
{
    return M(0,0);
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,2,2> &M)
{
    return M(0,0)*M(1,1) - M(1,0)*M(0,1);
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,3,3> &M)
{
    return ( M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) )
         - ( M(2,0)*M(1,1)*M(0,2) + M(2,1)*M(1,2)*M(0,0) + M(2,2)*M(1,0)*M(0,1) );
}

template<> inline ArrayDbl determinant(FMT_UNUSED const Eigen::Matrix<ArrayDbl,0,3> &M)
{
    return ArrayDbl();
}

template<> inline ArrayDbl determinant(FMT_UNUSED const Eigen::Matrix<ArrayDbl,3,0> &M)
{
    return ArrayDbl();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,1,2> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,2,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,1,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,3,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,2,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline ArrayDbl determinant(const Eigen::Matrix<ArrayDbl,3,2> &M)
{
	return determinant(normal_matrix(M)).sqrt();
}


/**
 * @brief Calculates inverse of rectangular matrix or pseudoinverse of non-rectangular matrix.
 */
template<int m, int n>
Eigen::Matrix<ArrayDbl,n,m> inverse(const Eigen::Matrix<ArrayDbl,m,n> &A) {
    // only for cases m > n
    return inverse(normal_matrix(A)) * A.transpose();
}


template<> inline Eigen::Matrix<ArrayDbl,1,1> inverse<1,1>(const Eigen::Matrix<ArrayDbl,1,1> &A)
{
	Eigen::Matrix<ArrayDbl,1,1> B;
    B(0,0) = A(0,0).inverse(); // 1/A(0,0)
    return B;
}

template<> inline Eigen::Matrix<ArrayDbl,2,2> inverse<2,2>(const Eigen::Matrix<ArrayDbl,2,2> &A)
{
	Eigen::Matrix<ArrayDbl,2,2> B;
	ArrayDbl det = determinant(A);

    B(0,0) = A(1,1) / det;
    B(0,1) = -A(0,1) / det;
    B(1,0) = -A(1,0) / det;
    B(1,1) = A(0,0) / det;
    return B;
}

template<> inline Eigen::Matrix<ArrayDbl,3,3> inverse<3,3>(const Eigen::Matrix<ArrayDbl,3,3> &A)
{
    Eigen::Matrix<ArrayDbl,3,3> B;

    B(0,0) = A(1,1)*A(2,2) - A(2,1)*A(1,2);
    B(1,0) = A(2,0)*A(1,2) - A(1,0)*A(2,2);
    B(2,0) = A(1,0)*A(2,1) - A(2,0)*A(1,1);

    ArrayDbl det = A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0);
    B(0,0) = B(0,0) / det;
    B(1,0) = B(1,0) / det;
    B(2,0) = B(2,0) / det;

    B(0,1) = (A(2,1)*A(0,2) - A(0,1)*A(2,2)) / det;
    B(1,1) = (A(0,0)*A(2,2) - A(2,0)*A(0,2)) / det;
    B(2,1) = (A(2,0)*A(0,1) - A(0,0)*A(2,1)) / det;

    B(0,2) = (A(0,1)*A(1,2) - A(1,1)*A(0,2)) / det;
    B(1,2) = (A(1,0)*A(0,2) - A(0,0)*A(1,2)) / det;
    B(2,2) = (A(0,0)*A(1,1) - A(1,0)*A(0,1)) / det;

    return B;
}

template<> inline Eigen::Matrix<ArrayDbl,2,1> inverse<1,2>(const Eigen::Matrix<ArrayDbl,1,2> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<ArrayDbl,3,1> inverse<1,3>(const Eigen::Matrix<ArrayDbl,1,3> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<ArrayDbl,3,2> inverse<2,3>(const Eigen::Matrix<ArrayDbl,2,3> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}


} // closing namespace eigen_tools

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


#endif /* PATCH_DATA_TABLE_HH_ */
