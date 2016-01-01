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
 * @file    bounding_box.cc
 * @brief   
 */

#include "system/system.hh"
#include "mesh/bounding_box.hh"


const double BoundingBox::epsilon = 64*numeric_limits<double>::epsilon();


BoundingBox::BoundingBox(const vector<Point> &points) {
	ASSERT_LESS( 0, points.size() );

	auto it = points.begin();
	max_vertex_ = min_vertex_ = *it;
	++it;
	for(; it != points.end(); ++it) expand( *it );
}

