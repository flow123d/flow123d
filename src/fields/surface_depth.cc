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
 * @file    surface_depth.cc
 * @brief
 */

#include "fields/surface_depth.hh"
#include "mesh/mesh.h"


SurfaceDepth::SurfaceDepth(const Mesh *mesh, std::string surface_region)
{
	bih_tree_ = new BIHTree();
	bih_tree_->add_boxes( const_cast<Mesh*>(mesh)->get_element_boxes(surface_region) );
	bih_tree_->construct();
}
