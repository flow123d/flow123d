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
 * @file    math_fce.h
 * @brief   
 */

#ifndef MATH_H
#define MATH_H

// TODO: Only mh_fe_values depends on this file. This file can be removed.

//! small matrix types
typedef double SmallVec1_t[1];
typedef double SmallVec2_t[2];
typedef double SmallVec3_t[3];
typedef double SmallVec4_t[4];

typedef SmallVec1_t *SmallMtx1;
typedef SmallVec2_t *SmallMtx2;
typedef SmallVec3_t *SmallMtx3;
typedef SmallVec4_t *SmallMtx4;


//! Numerical helpers
#include <float.h>
// DBL_MIN approx.= 2.22507e-308, least nonzero double

// this should be used when we want to keep first digit nonzero
// such number have full double precision, smaller numbers loose precision
#define NUM_ZERO DBL_MIN/DBL_EPSILON

// each use of ESP(value) indicate that we use some small value
// this small value should be relative to the number we want to compare
#define EPS(value) (value)

// OBSOLETE
#define ZERO		EPS(1e-12)

#define DBL_EQ(i,j) (fabs((i)-(j))<NUM_ZERO)
#define DBL_GE(i,j) ((i)>(j)-NUM_ZERO)
#define DBL_LE(i,j) ((i)<(j)+NUM_ZERO)
#define DBL_GT(i,j) ((i)>(j)+NUM_ZERO)
#define DBL_LT(i,j) ((i)<(j)-NUM_ZERO)

//! Usefull math macros
#define SQUARE(x)   ((x) * (x))
#define SGN(x)      ( ((x)>ZERO)? (1) :  ( ((x)<(-ZERO))? (-1) : (0) ) )
#define SUBDET2(i,j,k,l) (a[(i)][(k)]*a[(j)][(l)]-a[(i)][(l)]*a[(j)][(k)])
#ifndef M_PI
    #define M_PI 3.14159265358979323846264338327950288f
#endif

// vector functions
double vector_length(double[]);
double scalar_product(double[],double[]);
void normalize_vector(double[]);
void scale_vector(double[],double);

void vector_product(double[],double[],double[]);
void vector_difference(double[],double[],double[]);
// small matrix operations
double Det3( SmallMtx3 a);
double MatrixInverse(double *a,double *b,int size);
void PrintSmallMatrix(double *mtx, int size );
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
