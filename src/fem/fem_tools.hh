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
	if ( fabs(det) < 4*std::numeric_limits<double>::epsilon() ) return B;

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
	B(0,1) = -A(1,0)*A(2,2) + A(2,0)*A(1,2);
	B(0,2) = A(1,0)*A(2,1) - A(2,0)*A(1,1);

	double det = A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2);
	if ( fabs(det) < 4*std::numeric_limits<double>::epsilon() ) return B;
	B(0,0)/=det; B(0,1)/=det; B(0,2)/=det;

	B(1,0) = (-A(0,1)*A(2,2) + A(2,1)*A(0,2)) / det;
	B(1,1) = (A(0,0)*A(2,2) - A(2,0)*A(0,2)) / det;
	B(1,2) = (-A(0,0)*A(2,1) + A(2,0)*A(0,1)) / det;

	B(2,0) = (A(0,1)*A(1,2) - A(1,1)*A(0,2)) / det;
	B(2,1) = (-A(0,0)*A(1,2) + A(1,0)*A(0,2)) / det;
	B(2,2) = (A(0,0)*A(1,1) - A(1,0)*A(0,1)) / det;

	return B;
}


/**
 * @brief Calculates inverse matrix of a non-rectangular matrix.
 */
//template<arma::uword m, arma::uword n>
//arma::mat::fixed<n,m> pinverse(const arma::mat::fixed<m,n> &A) {
//    return A.t() * inverse(normal_matrix(A));
//}






#endif /* MAPPING_HH_ */
