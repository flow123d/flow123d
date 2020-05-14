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
 * @file    qmidpoint.hh
 * @brief   Midpoint rule qudrature.
 * @author  Pavel Exner
 */

#ifndef QMIDPOINT_HH_
#define QMIDPOINT_HH_

#include "system/global_defs.h"
#include "quadrature/quadrature.hh"
#include "mesh/point.hh"
#include "system/armor.hh"


/** @brief Class representing midpoint rule, with uniformly distributed points of the same weight.
 */
class QMidpoint : public Quadrature {
public:
    /// Empty constructor
    QMidpoint(unsigned int n_quadrature_points)
    : Quadrature(1, n_quadrature_points) {
        
        double qweight = 1.0/n_quadrature_points;
        for(unsigned int q=0; q < n_quadrature_points; q++) {
            this->weights[q] = qweight;
            this->set(q) = Armor::ArmaVec<double, 1>({(0.5+ q)*qweight});
        }
    }
};

#endif // QMIDPOINT_HH_
