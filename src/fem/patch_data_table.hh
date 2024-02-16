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
 * @file    patch_data_table.hh
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#ifndef PATCH_DATA_TABLE_HH_
#define PATCH_DATA_TABLE_HH_

#include <Eigen/Core>
#include <Eigen/Dense>


template <unsigned int N>
using ColData = Eigen::Array<double,N,1>;

namespace eigen_tools {

typedef ColData<300> Vec300;



/**
 * @brief Calculates determinant of a rectangular matrix.
 */
template<class T>
Vec300 determinant(const T &M);



inline Eigen::Matrix<Vec300,1,1> normal_matrix(const Eigen::Matrix<Vec300,1,2> &A) {
	Eigen::Matrix<Vec300,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1);
    return res;
}

inline Eigen::Matrix<Vec300,1,1> normal_matrix(const Eigen::Matrix<Vec300,2,1> &A) {
	Eigen::Matrix<Vec300,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0);
    return res;
}

inline Eigen::Matrix<Vec300,1,1> normal_matrix(const Eigen::Matrix<Vec300,1,3> &A) {
	Eigen::Matrix<Vec300,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    return res;
}

inline Eigen::Matrix<Vec300,1,1> normal_matrix(const Eigen::Matrix<Vec300,3,1> &A) {
	Eigen::Matrix<Vec300,1,1> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    return res;
}

inline Eigen::Matrix<Vec300,2,2> normal_matrix(const Eigen::Matrix<Vec300,2,3> &A) {
    Eigen::Matrix<Vec300,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(0,1)*A(0,1)+A(0,2)*A(0,2);
    res(0,1) = A(0,0)*A(1,0)+A(0,1)*A(1,1)+A(0,2)*A(1,2);
    res(1,0) = A(1,0)*A(0,0)+A(1,1)*A(0,1)+A(1,2)*A(0,2);
    res(1,1) = A(1,0)*A(1,0)+A(1,1)*A(1,1)+A(1,2)*A(1,2);
    return res;
}

inline Eigen::Matrix<Vec300,2,2> normal_matrix(const Eigen::Matrix<Vec300,3,2> &A) {
	Eigen::Matrix<Vec300,2,2> res;
    res(0,0) = A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    res(0,1) = A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    res(1,0) = A(0,1)*A(0,0)+A(1,1)*A(1,0)+A(2,1)*A(2,0);
    res(1,1) = A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    return res;
}



template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,1,1> &M)
{
    return M(0,0);
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,2,2> &M)
{
    return M(0,0)*M(1,1) - M(1,0)*M(0,1);
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,3,3> &M)
{
    return ( M(0,0)*M(1,1)*M(2,2) + M(0,1)*M(1,2)*M(2,0) + M(0,2)*M(1,0)*M(2,1) )
         - ( M(2,0)*M(1,1)*M(0,2) + M(2,1)*M(1,2)*M(0,0) + M(2,2)*M(1,0)*M(0,1) );
}

template<> inline Vec300 determinant(FMT_UNUSED const Eigen::Matrix<Vec300,0,3> &M)
{
    return Vec300();
}

template<> inline Vec300 determinant(FMT_UNUSED const Eigen::Matrix<Vec300,3,0> &M)
{
    return Vec300();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,1,2> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,2,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,1,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,3,1> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,2,3> &M)
{
    return determinant(normal_matrix(M)).sqrt();
}

template<> inline Vec300 determinant(const Eigen::Matrix<Vec300,3,2> &M)
{
	return determinant(normal_matrix(M)).sqrt();
}


/**
 * @brief Calculates inverse of rectangular matrix or pseudoinverse of non-rectangular matrix.
 */
template<int m, int n>
Eigen::Matrix<Vec300,n,m> inverse(const Eigen::Matrix<Vec300,m,n> &A) {
    // only for cases m > n
    return inverse(normal_matrix(A)) * A.transpose();
}


template<> inline Eigen::Matrix<Vec300,1,1> inverse<1,1>(const Eigen::Matrix<Vec300,1,1> &A)
{
	Eigen::Matrix<Vec300,1,1> B;
    B(0,0) = A(0,0).inverse(); // 1/A(0,0)
    return B;
}

template<> inline Eigen::Matrix<Vec300,2,2> inverse<2,2>(const Eigen::Matrix<Vec300,2,2> &A)
{
	Eigen::Matrix<Vec300,2,2> B;
	Vec300 det = determinant(A);

    B(0,0) = A(1,1) / det;
    B(0,1) = -A(0,1) / det;
    B(1,0) = -A(1,0) / det;
    B(1,1) = A(0,0) / det;
    return B;
}

template<> inline Eigen::Matrix<Vec300,3,3> inverse<3,3>(const Eigen::Matrix<Vec300,3,3> &A)
{
    Eigen::Matrix<Vec300,3,3> B;

    B(0,0) = A(1,1)*A(2,2) - A(2,1)*A(1,2);
    B(0,1) = A(2,0)*A(1,2) - A(1,0)*A(2,2);
    B(0,2) = A(1,0)*A(2,1) - A(2,0)*A(1,1);

    Vec300 det = A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
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

template<> inline Eigen::Matrix<Vec300,2,1> inverse<1,2>(const Eigen::Matrix<Vec300,1,2> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<Vec300,3,1> inverse<1,3>(const Eigen::Matrix<Vec300,1,3> &A)
{
    return A.transpose() * inverse(normal_matrix(A));
}

template<> inline Eigen::Matrix<Vec300,3,2> inverse<2,3>(const Eigen::Matrix<Vec300,2,3> &A)
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
