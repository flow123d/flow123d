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
 * @file    boundaries.cc
 * @ingroup mesh
 * @brief   Boundary conditions
 */

#include "system/system.hh"
#include "system/xio.h"
#include "mesh/mesh.h"
#include "mesh/boundaries.h"
#include "mesh/accessors.hh"


flow::VectorId<unsigned int> Boundary::id_to_bcd;


Boundary::Boundary()
: edge_idx_(Mesh::undef_idx), bc_ele_idx_(Mesh::undef_idx),
  mesh_(NULL)
{}


Element * Boundary::element() {
    return &( mesh_->bc_elements[bc_ele_idx_] );
}

Edge * Boundary::edge() {
    return &( mesh_->edges[edge_idx_] );
}

ElementAccessor<3> Boundary::element_accessor()
{
	return mesh_->element_accessor(bc_ele_idx_, true);
}

//-----------------------------------------------------------------------------
// vim: set cindent:

