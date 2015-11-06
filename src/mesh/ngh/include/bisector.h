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
 * @file    bisector.h
 * @brief   
 */

#ifndef bisectorH
#define bisectorH

#include <iostream>

#include "system/global_defs.h"
#include "mesh/element_impls.hh"

#include "point.h"
#include "myvector.h"

class TBisector {
protected:
    static int numberInstance;
    int id;

    TPoint* X0;
    TVector* U;

    int generateId();

public:
    TBisector();
    TBisector(const TPoint&, const TVector&);
    TBisector(const TPoint&, const TPoint&);
    TBisector(const TBisector &);
    TBisector(const Element &);
    ~TBisector();

    TBisector & operator =(const TBisector&);
    friend std::ostream & operator <<(std::ostream&, const TBisector&);

    void SetPoints(const TPoint&, const TPoint&);
    bool Belong(const TPoint&) const;

    void SetPoint(const TPoint&);
    const TPoint &GetPoint() const;
    TPoint GetPoint(double) const;
    void GetParameter(const TPoint&, double &, bool &) const;

    void SetVector(const TVector&);
    const TVector &GetVector() const;

    static int getNumInstances() {
        return TBisector::numberInstance;
    }
};

#endif
