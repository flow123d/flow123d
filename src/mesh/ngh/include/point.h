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
 * @file    point.h
 * @brief   
 */

#ifndef pointH
#define pointH

#include <iostream>

#include "myvector.h"

class TVector;

class TPoint {
private:
    static int numberInstance;
    int id;

    double x;
    double y;
    double z;

    int generateId();

public:
    TPoint();
    TPoint(double, double, double);
    TPoint(const TPoint&);
    ~TPoint();

    TPoint * operator =(TPoint*);
    TPoint * operator +(TPoint*);
    bool operator ==(TPoint*);

    TPoint & operator =(const TPoint&);
    TPoint & operator =(const TVector&);
    TVector operator -(const TPoint&) const;
    TPoint operator +(const TPoint&) const;
    bool operator ==(const TPoint&) const;
    friend std::ostream & operator <<(std::ostream&, const TPoint&);

    void SetCoord(double, double, double);
    void SetCoord(const TVector&);

    double X() const;
    double Y() const;
    double Z() const;

    double Get(int) const;

    static int getNumInstances() {
        return TPoint::numberInstance;
    }
};

#endif
