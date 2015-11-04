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
 * @file    intersection.h
 * @brief   
 */

#ifndef intersectionH
#define intersectionH

#include "bisector.h"
#include "abscissa.h"
#include "point.h"
#include "plain.h"
#include "triangle.h"
#include "tetrahedron.h"
#include "intersectionLocal.h"

typedef enum Intersections{
        none,
        unknown,
        point,
        line,
        area
} TIntersectionType;

typedef enum Positions{
        skew,
        parallel,
        intersecting,
        same,
        belong
} TPosition;

void GetIntersection(const TBisector &, const TBisector &, TPosition &,
                     double &, double &);
void GetIntersection(const TAbscissa &, const TAbscissa &, TPosition &,
                     double &, double &);
void GetIntersection(const TAbscissa &, const TAbscissa &, IntersectionLocal * & insec);

void GetIntersection(const TBisector &, const TAbscissa &, TPosition &,
                     double &, double &);
void GetIntersection(const TBisector &, const TAbscissa &, IntersectionLocal * & insec);

void GetIntersection(const TAbscissa &, const TBisector &, TPosition &,
                     double &, double &);
void GetIntersection(const TAbscissa &, const TBisector &, IntersectionLocal * & insec);

void GetIntersection(const TPlain &, const TPlain &,
                     TPosition &, TBisector *);
void GetIntersection(const TPlain &, const TBisector &,
                     TPosition &, TPoint *);
void GetIntersection(const TBisector &, const TPlain &,
                     TPosition &, double &);
void GetIntersection(const TBisector &, const TPlain &,
                     TPosition &, TPoint *);
void GetIntersection(const TTriangle &, const TTriangle &,
                     TIntersectionType &, double &);
void GetIntersection(const TBisector &, const TTriangle &, IntersectionLocal * & insec);

void GetIntersection(const TAbscissa &, const TTriangle &, IntersectionLocal * & insec);

void GetIntersection(const TAbscissa &, const TTetrahedron &,
                     TIntersectionType &, double &);
void GetIntersection(const TTriangle &, const TTetrahedron &,
                     TIntersectionType &, double &);

template<class A, class B> bool QuickIntersectionTest(const A &a, const B &b);

double Distance(const TBisector &, const TPoint &);
double Distance(const TPlain &, const TPoint &);
double Distance(const TPoint &, const TPoint &);
#endif
