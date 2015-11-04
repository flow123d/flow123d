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
 * @file    mathfce.h
 * @brief   
 */

#ifndef mainH
#define mainH

#include <limits>

const double epsilon = std::numeric_limits<double>::epsilon();


namespace mathfce{
        bool IsZero(double);
        bool IsEqual(double, double);
}
        double Determinant3(double [3][3]);
#endif
