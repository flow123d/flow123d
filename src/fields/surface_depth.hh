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
 * @file    surface_depth.hh
 * @brief
 */

#ifndef SURFACE_DEPTH_HH_
#define SURFACE_DEPTH_HH_

#include <string>
#include "mesh/bih_tree.hh"

class Mesh;


class SurfaceDepth
{
public:
	/**
	 * Constructor
	 */
	SurfaceDepth(const Mesh *mesh, std::string surface_region, std::string surface_direction);

	/**
	 * Destructor
	 */
	~SurfaceDepth();

protected:
	/// Tree of mesh elements
	BIHTree * bih_tree_;
};

#endif /* SURFACE_DEPTH_HH_ */
