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
 * @file    mathfce.cpp
 * @brief   
 */

#include <cmath>
#include "mesh/ngh/include/mathfce.h"


bool mathfce::IsZero(double x)
{
  if (fabs(x) < epsilon)
    return true;
  return false;
}

bool mathfce::IsEqual(double x, double y)
{
  if (fabs(x - y) < epsilon)
    return true;
  return false;
}


double Determinant3( double a[ 3 ][ 3 ] )
{
  double rc;
  rc =   a[ 0 ][ 0 ] * a[ 1 ][ 1 ] * a[ 2 ][ 2 ]
       + a[ 1 ][ 0 ] * a[ 2 ][ 1 ] * a[ 0 ][ 2 ]
       + a[ 2 ][ 0 ] * a[ 0 ][ 1 ] * a[ 1 ][ 2 ]
       - a[ 0 ][ 2 ] * a[ 1 ][ 1 ] * a[ 2 ][ 0 ]
       - a[ 1 ][ 2 ] * a[ 2 ][ 1 ] * a[ 0 ][ 0 ]
       - a[ 2 ][ 2 ] * a[ 0 ][ 1 ] * a[ 1 ][ 0 ];
  return rc;
}

