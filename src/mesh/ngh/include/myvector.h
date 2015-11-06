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
 * @file    myvector.h
 * @brief   
 */

#ifndef vectorH
#define vectorH

class TPoint;


class TVector {
protected:
    static int numberInstance;
    int id;

    double coors[ 3 ];
    double length;

    int generateId();

    void Compute();
    void CompLength();

public:
    TVector();
    TVector(double, double, double);
    TVector(TPoint, TPoint);
    TVector(const TVector &);
    ~TVector();

    double Length() const;
    void Get(double&, double&, double&) const;
    void Get(double*) const;
    double Get(int) const;
    void SetVector(double, double, double);
    bool IsZero();

    double X1() const;
    double X2() const;
    double X3() const;

    TVector & operator =(const TPoint&);
    TVector operator +(const TVector&);
    TVector operator +(const TPoint&);
    TVector operator -(const TVector&);
    friend TVector operator*(const TVector&, double);
    friend TVector operator*(double, const TVector&);
    bool operator ==(const TVector&);

    static int getNumInstances() {
        return TVector::numberInstance;
    }
};

TVector Cross(const TVector&, const TVector&);
double Dot(const TVector&, const TVector&);
bool AreParallel(const TVector&, const TVector&);
bool ArePerpendicular(const TVector&, const TVector&);

#endif
