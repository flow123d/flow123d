/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Boundary conditions
 * @ingroup mesh
 *
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

