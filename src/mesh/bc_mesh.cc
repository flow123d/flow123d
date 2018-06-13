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
 * @file    bc_mesh.cc
 * @ingroup mesh
 * @brief   Mesh construction
 */



#include "mesh/bc_mesh.hh"
#include "mesh/long_idx.hh"
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"
#include "mesh/neighbours.h"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"



BCMesh::BCMesh(Mesh* parent_mesh)
: parent_mesh_(parent_mesh)
{
	this->init_element_vector(0);
	this->init_node_vector(0);
}


BCMesh::~BCMesh()
{}


Range<ElementAccessor<3>> BCMesh::elements_range(bool boundary) const
{
    return Range<ElementAccessor<3>>(parent_mesh_, parent_mesh_->bulk_size_, parent_mesh_->element_vec_.size());
}


unsigned int BCMesh::n_elements(bool boundary) const {
	return parent_mesh_->element_ids_.size()-parent_mesh_->bulk_size_;
}


unsigned int BCMesh::elements_shift(bool boundary) const {
	return parent_mesh_->bulk_size_;
}


unsigned int BCMesh::n_edges() const
{ return parent_mesh_->edges.size(); }


unsigned int BCMesh::n_vb_neighbours() const
{ return parent_mesh_->vb_neighbours_.size(); }


Edge *BCMesh::edge(unsigned int i) {
    return parent_mesh_->edge(i);
}


Neighbour *BCMesh::vb_neighbour(unsigned int i) {
    return parent_mesh_->vb_neighbour(i);
}


Partitioning *BCMesh::get_part()
{ return parent_mesh_->get_part(); }


Distribution *BCMesh::get_el_ds() const
{ return parent_mesh_->el_ds; }


LongIdx *BCMesh::get_row_4_el() const
{ return parent_mesh_->row_4_el; }


LongIdx *BCMesh::get_el_4_loc() const
{ return parent_mesh_->el_4_loc; }

