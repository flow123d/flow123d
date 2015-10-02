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
 * @file    tetrahedron.h
 * @brief   
 */

#ifndef tetrahedronH
#define tetrahedronH

#include "point.h"
#include "triangle.h"

class TTetrahedron {
private:
    static int numberInstance;
    int id;

    TPoint X1;
    TPoint X2;
    TPoint X3;
    TPoint X4;

    TTriangle T1;
    TTriangle T2;
    TTriangle T3;
    TTriangle T4;

    TAbscissa* A1;
    TAbscissa* A2;
    TAbscissa* A3;
    TAbscissa* A4;
    TAbscissa* A5;
    TAbscissa* A6;

    double volume;

    int generateId();
    void ComputeVolume();

public:
    TTetrahedron();
    TTetrahedron(const TPoint&, const TPoint&, const TPoint&, const TPoint&);
    ~TTetrahedron();

    const TTriangle &GetTriangle(int) const;
    const TAbscissa &GetAbscissa(int) const;
    const TPoint &GetPoint(int) const;

    double GetMin(int) const;
    double GetMax(int) const;

    double GetVolume();

    void SetPoints(const TPoint& P1, const TPoint& P2, const TPoint& P3, const TPoint& P4);
    bool IsInner(const TPoint&) const;

    static int getNumInstances() {
        return TTetrahedron::numberInstance;
    }
};
#endif
