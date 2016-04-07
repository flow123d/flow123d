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
 * @file    bisector.cpp
 * @brief   
 */

#include <iostream>

#include "system/global_defs.h"
#include "mesh/elements.h"

#include "mesh/ngh/include/bisector.h"
#include "mesh/ngh/include/mathfce.h"
#include "mesh/ngh/include/intersection.h"

using namespace mathfce;

int TBisector::numberInstance = 0;

int TBisector::generateId() {
    return TBisector::numberInstance++;
}

TBisector::TBisector() {
    id = generateId();

    X0 = new TPoint(0, 0, 0);
    U = new TVector(0, 0, 0);
}

TBisector::TBisector(const TPoint &XX0, const TVector &UU) {
    id = generateId();

    X0 = new TPoint(XX0);
    U = new TVector(UU);
}

TBisector::TBisector(const TPoint &P0, const TPoint &P1) {
    id = generateId();

    X0 = new TPoint(P0);
    U = new TVector(P0, P1);
}

TBisector::TBisector(const Element & ele) {
    id = generateId();
    OLD_ASSERT_EQUAL(ele.dim(), 1);

    X0 = new TPoint(ele.node[0]->point()(0), ele.node[0]->point()(1), ele.node[0]->point()(2));
    U = new TVector(*X0, TPoint(ele.node[1]->point()(0), ele.node[1]->point()(1), ele.node[1]->point()(2)) );
}

TBisector::TBisector(const  TBisector &x)
{
    id = generateId();

    X0 = new TPoint(*(x.X0));
    U = new TVector(*(x.U));
}

TBisector::~TBisector() {
    delete X0;
    delete U;
}

TBisector & TBisector::operator =(const TBisector &b) {
    *(*this).U = *b.U;
    *(*this).X0 = *b.X0;

    return *this;
}

std::ostream & operator <<(std::ostream &stream, const TBisector &b) {
    stream << "U = (" << b.U->X1() << ", " << b.U->X2() << ", " << b.U->X3() << ")\n";
    stream << "X0 = [" << b.X0->X() << ", " << b.X0->Y() << ", " << b.X0->Z() << "]\n";

    return stream;
}

void TBisector::SetPoint(const TPoint &P) {
    *X0 = P;
}

void TBisector::SetVector(const TVector &UU) {
    *U = UU;
}

void TBisector::SetPoints(const TPoint &P0, const TPoint &P1) {
    *X0 = P0;
    *U = (TPoint) P1 - P0;
}

const TVector &TBisector::GetVector() const {

    return *(this->U);
}

const TPoint &TBisector::GetPoint() const
{
    return *(this->X0);
}

TPoint TBisector::GetPoint(double t) const {
    TPoint tmp;

    tmp.SetCoord(t * *U);
    tmp = tmp + *X0;

    return tmp;
}

void TBisector::GetParameter(const TPoint &P, double &t, bool &onBisector) const {
	t = (P.X() - X0->X()) / U->X1();
	onBisector = (fabs( (P.Y() - X0->Y()) / U->X2() - t ) < epsilon) & (fabs( (P.Z() - X0->Z()) / U->X3() - t ) < epsilon);
}

bool TBisector::Belong(const TPoint &P) const {
    if (IsZero(Distance(*this, P))) {
        return true;
    } else {
        return false;
    }
}
