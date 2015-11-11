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
 * @file    triangle.h
 * @brief   
 */

#ifndef triangleH
#define triangleH

#include "point.h"
#include "plain.h"
#include "abscissa.h"
#include "mesh/bounding_box.hh"
#include "mesh/elements.h"

class TTriangle {
private:
    static int numberInstance;
    int id;

    TPoint X1;
    TPoint X2;
    TPoint X3;

    TAbscissa* A1;
    TAbscissa* A2;
    TAbscissa* A3;

    TPlain* pl;

    BoundingBox boundingBox;

    double area;

    int generateId();

    void ComputeArea();

public:
    TTriangle();
    TTriangle(const TTriangle&);
    TTriangle(const TPoint&, const TPoint&, const TPoint&);
    TTriangle(const Element&);
    void update();

    ~TTriangle();

    TTriangle & operator =(const TTriangle &t);

    const TPlain &GetPlain() const;
    const TAbscissa &GetAbscissa(int) const;
    const TPoint &GetPoint(int) const;

    void SetPoints(const TPoint&, const TPoint&, const TPoint&);

    double GetMin(int) const;
    double GetMax(int) const;

    double GetArea();
    BoundingBox &get_bounding_box();

    bool IsInner(const TPoint&) const;

    static int getNumInstances() {
        return TTriangle::numberInstance;
    }
};

#endif

