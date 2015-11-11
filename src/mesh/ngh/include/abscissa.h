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
 * @file    abscissa.h
 * @brief   
 */

#ifndef abscissaH
#define abscissaH

#include "bisector.h"
#include "point.h"
#include "mesh/bounding_box.hh"

class TAbscissa : public TBisector {
private:
    static int numberInstance;
    int id;

    BoundingBox boundingBox;

    double length;

    int generateId();
    void ComputeLength();

public:
    TAbscissa();
    TAbscissa(const TPoint&, const TPoint&);
    TAbscissa(const Element&);
    ~TAbscissa();

    TAbscissa & operator =(const TAbscissa&);

    double Length();
    BoundingBox &get_bounding_box();

    void SetPoints(const TPoint&, const TPoint&);

    double GetMin(int) const;
    double GetMax(int) const;

    static int getNumInstances() {
        return TAbscissa::numberInstance;
    }
};

#endif
