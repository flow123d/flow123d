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
 * @file    plain.cpp
 * @brief   
 */

#include "mesh/ngh/include/plain.h"
#include "mesh/ngh/include/mathfce.h"
#include "mesh/ngh/include/intersection.h"

using namespace mathfce;

int TPlain::numberInstance = 0;

int TPlain::generateId() {
    return TPlain::numberInstance++;
}

TPlain::TPlain() {
    id = generateId();

    U = new TVector(0, 0, 0);
    V = new TVector(0, 0, 0);
    N = new TVector(0, 0, 0);
    X = new TPoint(0, 0, 0);

    Compute();
}

TPlain::TPlain(const TVector &UU, const TVector &VV, const TPoint &XX) {
    id = generateId();

    U = new TVector(UU);
    V = new TVector(VV);
    N = new TVector(0, 0, 0);
    X = new TPoint(XX);

    Compute();
}

TPlain::TPlain(const TPoint &P1, const TPoint &P2, const TPoint &P3) {
    id = generateId();

    U = new TVector(P1, P2);
    V = new TVector(P1, P3);
    N = new TVector(0, 0, 0);
    X = new TPoint(P1);

    Compute();
}

TPlain::TPlain(const TPlain &p) {
    id = generateId();

    N = new TVector();
    U = new TVector();
    V = new TVector();
    X = new TPoint();

    a = p.a;
    b = p.b;
    c = p.c;
    d = p.d;

    *N = *p.N;
    *U = *p.U;
    *V = *p.V;
    *X = *p.X;
}

void TPlain::Compute() {
    *N = Cross(*U, *V);
    a = N->X1();
    b = N->X2();
    c = N->X3();
    d = -a * X->X() - b * X->Y() - c * X->Z();
    return;
}

TPlain::~TPlain() {
    delete U;
    delete V;
    delete N;
    delete X;
}

double TPlain::GetA() const {
    return a;
}

double TPlain::GetB() const {
    return b;
}

double TPlain::GetC() const {
    return c;
}

double TPlain::GetD() const {
    return d;
}

const TVector &TPlain::GetNormal() const {
    return *N;
}

const TVector &TPlain::GetU() const {
    return *U;
}

const TVector &TPlain::GetV() const {
    return *V;
}

const TPoint &TPlain::GetPoint() const {
    return *X;
}

TPoint TPlain::GetPoint(double r, double s) const {
    TPoint tmp;
    tmp = r * *U + s * *V + *X;
    return tmp;
}

bool TPlain::Belong(const TPoint &P) const {
    if (IsZero(Distance(*this, P)))
        return true;
    return false;
}

TPlain & TPlain::operator =(const TPlain &p) {
    a = p.a;
    b = p.b;
    c = p.c;
    d = p.d;
    *(*this).U = *p.U;
    *(*this).V = *p.V;
    *(*this).N = *p.N;
    *(*this).X = *p.X;
    return *this;
}

void TPlain::SetPoints(const TPoint &P1, const TPoint &P2, const TPoint &P3) {
    *U = (TPoint) P2 - P1;
    *V = (TPoint) P3 - P1;
    *X = P1;
    Compute();
}

