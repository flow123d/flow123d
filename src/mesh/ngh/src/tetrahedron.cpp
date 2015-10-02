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
 * @file    tetrahedron.cpp
 * @brief   
 */

#include <cmath>
#include "system/exc_common.hh"

#include "mesh/ngh/include/tetrahedron.h"
#include "mesh/ngh/include/intersection.h"
#include "mesh/ngh/include/mathfce.h"

int TTetrahedron::numberInstance = 0;

int TTetrahedron::generateId() {
    return TTetrahedron::numberInstance++;
}

TTetrahedron::TTetrahedron() {
	id = generateId();

	A1 = new TAbscissa();
	A2 = new TAbscissa();
	A3 = new TAbscissa();
	A4 = new TAbscissa();
	A5 = new TAbscissa();
	A6 = new TAbscissa();

	volume = 0.0;
}

TTetrahedron::TTetrahedron(const TPoint& X1, const TPoint& X2,
        const TPoint& X3, const TPoint& X4) {
    id = generateId();

    this->X1 = X1;
    this->X2 = X2;
    this->X3 = X3;
    this->X4 = X4;

    T1.SetPoints(X2, X3, X4);
    T2.SetPoints(X1, X3, X4);
    T3.SetPoints(X1, X2, X4);
    T4.SetPoints(X1, X2, X3);

    A1 = new TAbscissa(X1, X2);
    A2 = new TAbscissa(X2, X3);
    A3 = new TAbscissa(X3, X1);
    A4 = new TAbscissa(X1, X4);
    A5 = new TAbscissa(X2, X4);
    A6 = new TAbscissa(X3, X4);

    ComputeVolume();
}

TTetrahedron::~TTetrahedron() {
    delete A1;
    delete A2;
    delete A3;
    delete A4;
    delete A5;
    delete A6;
}

const TTriangle &TTetrahedron::GetTriangle(int i) const {
    switch (i) {
        case 1: return T1;
            break;
        case 2: return T2;
            break;
        case 3: return T3;
            break;
        case 4: return T4;
            break;
        default: THROW( ExcAssertMsg() << EI_Message("Unknown number of the triangle of the tetrahedron.") );
    }
}

const TAbscissa &TTetrahedron::GetAbscissa(int i) const {
    switch (i) {
        case 1: return *A1;
            break;
        case 2: return *A2;
            break;
        case 3: return *A3;
            break;
        case 4: return *A4;
            break;
        case 5: return *A5;
            break;
        case 6: return *A6;
            break;
        default: THROW( ExcAssertMsg() << EI_Message("Unknown number of the triangle of the tetrahedron.") );
    }
}

const TPoint &TTetrahedron::GetPoint(int i) const {
    switch (i) {
        case 1: return X1;
            break;
        case 2: return X2;
            break;
        case 3: return X3;
            break;
        case 4: return X4;
            break;
        default: THROW( ExcAssertMsg() << EI_Message("Unknown number of the point of the tetrahedron.") );
    }
}

double TTetrahedron::GetMin(int i) const {
    double min = X1.Get(i);

    if (X2.Get(i) < min) {
        min = X2.Get(i);
    }
    if (X3.Get(i) < min) {
        min = X3.Get(i);
    }
    if (X4.Get(i) < min) {
        min = X4.Get(i);
    }

    return min;
}

double TTetrahedron::GetMax(int i) const {
    double max = X1.Get(i);

    if (X2.Get(i) > max) {
        max = X2.Get(i);
    }
    if (X3.Get(i) > max) {
        max = X3.Get(i);
    }
    if (X4.Get(i) > max) {
        max = X4.Get(i);
    }

    return max;
}

double TTetrahedron::GetVolume() {
    return volume;
}

void TTetrahedron::ComputeVolume() {
    double a[ 3 ][ 3 ];

    a[ 0 ][ 0 ] = X2.X() - X1.X();
    a[ 0 ][ 1 ] = X2.Y() - X1.Y();
    a[ 0 ][ 2 ] = X2.Z() - X1.Z();
    a[ 1 ][ 0 ] = X3.X() - X1.X();
    a[ 1 ][ 1 ] = X3.Y() - X1.Y();
    a[ 1 ][ 2 ] = X3.Z() - X1.Z();
    a[ 2 ][ 0 ] = X4.X() - X1.X();
    a[ 2 ][ 1 ] = X4.Y() - X1.Y();
    a[ 2 ][ 2 ] = X4.Z() - X1.Z();

    volume = fabs(Determinant3(a)) / 6.0;
}

void TTetrahedron::SetPoints(const TPoint& P1, const TPoint& P2, const TPoint& P3, const TPoint& P4) {
	X1 = P1;
	X2 = P2;
	X3 = P3;
	X4 = P4;

	T1.SetPoints(P2, P3, P4);
	T2.SetPoints(P1, P3, P4);
	T3.SetPoints(P1, P2, P4);
	T4.SetPoints(P1, P2, P3);

	A1->SetPoints(P1, P2);
	A2->SetPoints(P2, P3);
	A3->SetPoints(P3, P1);
	A4->SetPoints(P1, P4);
	A5->SetPoints(P2, P4);
	A6->SetPoints(P3, P4);

	ComputeVolume();
}

bool TTetrahedron::IsInner(const TPoint& P) const {
	TVector N, U1(X1, X2), U2(X1, X3), U3(X1, X4), U4(X2, X3), U5(X2, X4), U6(X2, X1);
	TVector Up1(X1, P), Up2(X2, P);

	N = Cross(U1, U2); //X4 is on opposite side of plain X1,X2,X3 than P
	if (Dot(N, U3) * Dot(N, Up1) < 0) {
		return false;
	}

	N = Cross(U1, U3); //X3 x P
	if (Dot(N, U2) * Dot(N, Up1) < 0) {
		return false;
	}

	N = Cross(U2, U3); //X2 x P
	if (Dot(N, U1) * Dot(N, Up1) < 0) {
		return false;
	}

	N = Cross(U4, U5); //X1 x P
	if (Dot(N, U6) * Dot(N, Up2) < 0) {
		return false;
	}

	return true;
}
