/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: quad.h,v 1.46 2011/04/15 08:10:07 zlb Exp $ */

/*
 * The parts of this file which are unnecessary for Flow123d have been commented out.
 */

#ifndef PHG_QUAD_H

typedef struct QUAD_ {
    const char *name;			/* name of the quadrature formulae */
    int dim;			/* dimension, 1: edge, 2: face, 3: tetra */
    int order;			/* exact for polynomials of order 'order' */
    int npoints;		/* number of points */
    FLOAT *points;		/* barycentric coordinates of quad. points */
    FLOAT *weights;		/* weights */
    SHORT id;			/* id (for use with reference count) */
} QUAD;

#define QUAD_DEFAULT	-1

/* 1D quadrature rules */
extern QUAD QUAD_1D_P1_;
#define QUAD_1D_P1 (&QUAD_1D_P1_)
extern QUAD QUAD_1D_P3_;
#define QUAD_1D_P2 (&QUAD_1D_P3_)
#define QUAD_1D_P3 (&QUAD_1D_P3_)
extern QUAD QUAD_1D_P5_;
#define QUAD_1D_P4 (&QUAD_1D_P5_)
#define QUAD_1D_P5 (&QUAD_1D_P5_)
extern QUAD QUAD_1D_P7_;
#define QUAD_1D_P6 (&QUAD_1D_P7_)
#define QUAD_1D_P7 (&QUAD_1D_P7_)
extern QUAD QUAD_1D_P9_;
#define QUAD_1D_P8 (&QUAD_1D_P9_)
#define QUAD_1D_P9 (&QUAD_1D_P9_)
extern QUAD QUAD_1D_P11_;
#define QUAD_1D_P10 (&QUAD_1D_P11_)
#define QUAD_1D_P11 (&QUAD_1D_P11_)
extern QUAD QUAD_1D_P13_;
#define QUAD_1D_P12 (&QUAD_1D_P13_)
#define QUAD_1D_P13 (&QUAD_1D_P13_)
extern QUAD QUAD_1D_P15_;
#define QUAD_1D_P14 (&QUAD_1D_P15_)
#define QUAD_1D_P15 (&QUAD_1D_P15_)
extern QUAD QUAD_1D_P17_;
#define QUAD_1D_P16 (&QUAD_1D_P17_)
#define QUAD_1D_P17 (&QUAD_1D_P17_)
extern QUAD QUAD_1D_P19_;
#define QUAD_1D_P18 (&QUAD_1D_P19_)
#define QUAD_1D_P19 (&QUAD_1D_P19_)
extern QUAD QUAD_1D_P21_;
#define QUAD_1D_P20 (&QUAD_1D_P21_)
#define QUAD_1D_P21 (&QUAD_1D_P21_)

/* 2D quadrature rules */
extern QUAD QUAD_2D_P1_;
#define QUAD_2D_P1 (&QUAD_2D_P1_)
extern QUAD QUAD_2D_P2_;
#define QUAD_2D_P2 (&QUAD_2D_P2_)
extern QUAD QUAD_2D_P3_;
#define QUAD_2D_P3 (&QUAD_2D_P3_)
extern QUAD QUAD_2D_P4_;
#define QUAD_2D_P4 (&QUAD_2D_P4_)
extern QUAD QUAD_2D_P5_;
#define QUAD_2D_P5 (&QUAD_2D_P5_)
extern QUAD QUAD_2D_P6_;
#define QUAD_2D_P6 (&QUAD_2D_P6_)
extern QUAD QUAD_2D_P7_;
#define QUAD_2D_P7 (&QUAD_2D_P7_)
extern QUAD QUAD_2D_P8_;
#define QUAD_2D_P8 (&QUAD_2D_P8_)
extern QUAD QUAD_2D_P9_;
#define QUAD_2D_P9 (&QUAD_2D_P9_)
extern QUAD QUAD_2D_P10_;
#define QUAD_2D_P10 (&QUAD_2D_P10_)
extern QUAD QUAD_2D_P11_;
#define QUAD_2D_P11 (&QUAD_2D_P11_)
extern QUAD QUAD_2D_P12_;
#define QUAD_2D_P12 (&QUAD_2D_P12_)
extern QUAD QUAD_2D_P13_;
#define QUAD_2D_P13 (&QUAD_2D_P13_)
extern QUAD QUAD_2D_P14_;
#define QUAD_2D_P14 (&QUAD_2D_P14_)
extern QUAD QUAD_2D_P15_;
#define QUAD_2D_P15 (&QUAD_2D_P15_)
extern QUAD QUAD_2D_P16_;
#define QUAD_2D_P16 (&QUAD_2D_P16_)
extern QUAD QUAD_2D_P17_;
#define QUAD_2D_P17 (&QUAD_2D_P17_)
extern QUAD QUAD_2D_P18_;
#define QUAD_2D_P18 (&QUAD_2D_P18_)
extern QUAD QUAD_2D_P19_;
#define QUAD_2D_P19 (&QUAD_2D_P19_)
extern QUAD QUAD_2D_P20_;
#define QUAD_2D_P20 (&QUAD_2D_P20_)
extern QUAD QUAD_2D_P21_;
#define QUAD_2D_P21 (&QUAD_2D_P21_)

/* 3D quadrature rules */
extern QUAD QUAD_3D_P1_;
#define QUAD_3D_P1 (&QUAD_3D_P1_)
extern QUAD QUAD_3D_P2_;
#define QUAD_3D_P2 (&QUAD_3D_P2_)
extern QUAD QUAD_3D_P3_;
#define QUAD_3D_P3 (&QUAD_3D_P3_)
extern QUAD QUAD_3D_P4_;
#define QUAD_3D_P4 (&QUAD_3D_P4_)
extern QUAD QUAD_3D_P5_;
#define QUAD_3D_P5 (&QUAD_3D_P5_)
extern QUAD QUAD_3D_P6_;
#define QUAD_3D_P6 (&QUAD_3D_P6_)
extern QUAD QUAD_3D_P7_;
#define QUAD_3D_P7 (&QUAD_3D_P7_)
extern QUAD QUAD_3D_P8_;
#define QUAD_3D_P8 (&QUAD_3D_P8_)
extern QUAD QUAD_3D_P9_;
#define QUAD_3D_P9 (&QUAD_3D_P9_)
extern QUAD QUAD_3D_P10_;
#define QUAD_3D_P10 (&QUAD_3D_P10_)
extern QUAD QUAD_3D_P11_;
#define QUAD_3D_P11 (&QUAD_3D_P11_)
extern QUAD QUAD_3D_P12_;
#define QUAD_3D_P12 (&QUAD_3D_P12_)
extern QUAD QUAD_3D_P13_;
#define QUAD_3D_P13 (&QUAD_3D_P13_)
extern QUAD QUAD_3D_P14_;
#define QUAD_3D_P14 (&QUAD_3D_P14_)

/*-------------------------- Permutation macros ----------------------------*/

/* 1D */
#define Perm2(a)	_F(a),_F(a)
#define Dup2(w)		_F(w)
#define Perm11(a)	_F(a),_F(1.)-(_F(a)), _F(1.)-(_F(a)),_F(a)
#define Dup11(w)	_F(w),_F(w)

/* 2D */
#define Perm3(a)	_F(a),_F(a),_F(a)
#define Dup3(w)		_F(w)
#define Perm21(a)	_F(a),_F(a),_F(1.)-(_F(a))-(_F(a)), \
			_F(a),_F(1.)-(_F(a))-(_F(a)),_F(a), \
			_F(1.)-(_F(a))-(_F(a)),_F(a),_F(a)
#define Dup21(w)	_F(w),_F(w),_F(w)
#define Perm111(a,b)	_F(a),_F(b),_F(1.)-(_F(a))-(_F(b)), \
			_F(a),_F(1.)-(_F(a))-(_F(b)),_F(b), \
			_F(b),_F(a),_F(1.)-(_F(a))-(_F(b)), \
			_F(b),_F(1.)-(_F(a))-(_F(b)),_F(a), \
			_F(1.)-(_F(a))-(_F(b)),_F(a),_F(b), \
			_F(1.)-(_F(a))-(_F(b)),_F(b),_F(a)
#define Dup111(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)

/* 3D */
#define Perm4(a)	_F(a),_F(a),_F(a),_F(a)
#define Dup4(w)		_F(w)
#define Perm31(a)	_F(a),_F(a),_F(a),_F(1.)-_F(3.)*(_F(a)), \
			_F(a),_F(a),_F(1.)-_F(3.)*(_F(a)),_F(a), \
			_F(a),_F(1.)-_F(3.)*(_F(a)),_F(a),_F(a), \
			_F(1.)-_F(3.)*(_F(a)),_F(a),_F(a),_F(a)
#define Dup31(w)	_F(w),_F(w),_F(w),_F(w)
#define Perm22(a)	_F(a),_F(a),_F(.5)-(_F(a)),_F(.5)-(_F(a)), \
			_F(a),_F(.5)-(_F(a)),_F(a),_F(.5)-(_F(a)), \
			_F(a),_F(.5)-(_F(a)),_F(.5)-(_F(a)),_F(a), \
			_F(.5)-(_F(a)),_F(a),_F(.5)-(_F(a)),_F(a), \
			_F(.5)-(_F(a)),_F(a),_F(a),_F(.5)-(_F(a)), \
			_F(.5)-(_F(a)),_F(.5)-(_F(a)),_F(a),_F(a)
#define Dup22(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm211(a,b)	_F(a),_F(a),_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(a),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b), \
			_F(a),_F(b),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(a),_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a), \
			_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(b), \
			_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b),_F(a), \
			_F(b),_F(a),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)), \
			_F(b),_F(a),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a), \
			_F(b),_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(a), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(a),_F(b), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(a),_F(b),_F(a), \
			_F(1.)-(_F(a))-(_F(a))-(_F(b)),_F(b),_F(a),_F(a)
#define Dup211(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w),\
			_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm0111(p,a,b,c) p,a,b,c, p,a,c,b, p,b,a,c, p,b,c,a, p,c,a,b, p,c,b,a
#define Perm1111(a,b,c) \
	Perm0111(_F(a),_F(b),_F(c),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(b),_F(a),_F(c),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(c),_F(a),_F(b),_F(1.)-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(1.)-(_F(a))-(_F(b))-(_F(c)),_F(a),_F(b),_F(c))
#define Dup1111(w)	Dup111(w), Dup111(w), Dup111(w), Dup111(w)


#define PHG_QUAD_H
#endif
