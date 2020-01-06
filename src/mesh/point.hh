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

#include "system/armor.hh"


/*
 * TODO:
 * - This Point notation is used in:
 * 		- fields, in particular value and value_list method
 * 		- io, output_element, output_mesh
 * 		- bih_tree, bounding box
 * 		- With own definition in fem/maps... !!
 * - not used in:
 *      - Mesh, Node
 *      -
 * Used single Point definition in mesh and reuse it elsewhere.
 * Do not use it in BIH, in order to make BIH separate library (WIP).
 * Try to use Armor instead of Armadillo fixed vectors.
 */

template <uint spacedim>
class Space {
public:
    typedef typename Armor::ArmaVec<double, spacedim> Point;
};


#endif /* POINT_HH_ */
