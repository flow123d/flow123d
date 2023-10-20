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






#endif /* MAPPING_HH_ */
