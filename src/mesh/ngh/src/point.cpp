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
 * @file    point.cpp
 * @brief   
 */

#include <iostream>

#include "system/exc_common.hh"
#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/mathfce.h"

using namespace mathfce;

int TPoint::numberInstance = 0;

int TPoint::generateId() {
    return TPoint::numberInstance++;
}

TPoint::TPoint() {
    x = 0;
    y = 0;
    z = 0;

    id = generateId();
}

TPoint::TPoint(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;

    id = generateId();
}

TPoint::TPoint(const TPoint& point) {
    x = point.X();
    y = point.Y();
    z = point.Z();

    id = generateId();
}

TPoint::~TPoint() {
}

TPoint& TPoint::operator =(const TPoint& P) {
    x = P.x;
    y = P.y;
    z = P.z;

    return *this;
}

bool TPoint::operator ==(const TPoint& P) const {
    if (IsEqual(x, P.x) && IsEqual(y, P.y) && IsEqual(z, P.z)) {
        return true;
    }
    return false;
}

TPoint& TPoint::operator =(const TVector& U) {
    x = U.X1();
    y = U.X2();
    z = U.X3();

    return *this;
}

TVector TPoint::operator -(const TPoint& P) const {
    return TVector( x - P.x, y - P.y, z - P.z);
}

TPoint TPoint::operator +(const TPoint& P) const {
    TPoint res;

    res.x = x + P.x;
    res.y = y + P.y;
    res.z = z + P.z;

    return res;
}

void TPoint::SetCoord(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

void TPoint::SetCoord(const TVector& U) {
    x = U.X1();
    y = U.X2();
    z = U.X3();
}

double TPoint::X() const {
    return x;
}

double TPoint::Y() const {
    return y;
}

double TPoint::Z() const {
    return z;
}

double TPoint::Get(int i) const {
    if (!(i >= 1 && i <= 3)) {
        THROW( ExcAssertMsg() << EI_Message( "Invalid specification of the element of the vector.") );
    }

    switch (i) {
        case 1: return x;
        case 2: return y;
        case 3: return z;
    }

    return 0.0;
}

std::ostream & operator <<(std::ostream& stream, const TPoint& P) {
    stream << "[" << P.x << " " << P.y << " " << P.z << "]";
    return stream;
}
