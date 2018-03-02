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

#include <armadillo>

#include "fields/surface_depth.hh"
#include "mesh/mesh.h"


SurfaceDepth::SurfaceDepth(const Mesh *mesh, std::string surface_region, std::string surface_direction)
{
	arma::vec3 surface_vec(surface_direction);
	ASSERT( surface_vec(2)!=0 ).error("Vertical plane of surface is forbidden!"); // may be not false
	surface_vec /= arma::norm(surface_vec, 2); // normalize to size == 1

	arma::vec3 ex, ey;
	ex(0) = 1; ex(1) = 0; ex(2) = - surface_vec(0) / surface_vec(2);
	ey(0) = 0; ey(1) = 1; ey(2) = - surface_vec(1) / surface_vec(2);

	arma::vec3 ta = ex - arma::dot(ex, surface_vec) * surface_vec;
	ta /= arma::norm(ta, 2); // normalize
	arma::vec3 tb = ey - arma::dot(ey, surface_vec) * surface_vec;
	tb /= arma::norm(tb, 2); // normalize

	bih_tree_ = new BIHTree();
	//bih_tree_->add_boxes( const_cast<Mesh*>(mesh)->get_element_boxes(surface_region) );
	//bih_tree_->construct();
}


SurfaceDepth::~SurfaceDepth()
{
	delete bih_tree_;
}
