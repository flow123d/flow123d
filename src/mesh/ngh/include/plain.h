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
 * @file    plain.h
 * @brief   
 */

#ifndef plainH
#define plainH

#include "myvector.h"
#include "point.h"

class TPlain {
private:
    static int numberInstance;
    int id;

    TVector* U;
    TVector* V;
    TVector* N;

    double a;
    double b;
    double c;
    double d;

    TPoint* X;

    int generateId();

    void Compute();

public:
    TPlain();
    TPlain(const TVector&, const TVector&, const TPoint&);
    TPlain(const TPoint&, const TPoint&, const TPoint&);
    TPlain(const TPlain&);
    ~TPlain();

    const TPoint &GetPoint() const;
    TPoint GetPoint(double, double) const;
    const TVector &GetNormal() const;
    const TVector &GetU() const;
    const TVector &GetV() const;
    double GetA() const;
    double GetB() const;
    double GetC() const;
    double GetD() const;
    bool Belong(const TPoint&) const;
    void SetPoints(const TPoint&, const TPoint&, const TPoint&);

    TPlain & operator =(const TPlain&);

    static int getNumInstances() {
        return TPlain::numberInstance;
    }
};

#endif
