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
 * @file    point.hh
 * @brief   
 */

#ifndef POINT_HH_
#define POINT_HH_

#include <armadillo>


/*
 * TODO:
 * need better resolution of various C++11 functionalities
 * e.g. following is supported from GCC 4.7
 */

template<int spacedim>
class Space {
public:
    typedef typename arma::vec::fixed<spacedim> Point;
};


#endif /* POINT_HH_ */
