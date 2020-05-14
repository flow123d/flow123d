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


#include "system/index_types.hh"
#include "mesh/bc_mesh.hh"
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"
#include "mesh/neighbours.h"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"



BCMesh::BCMesh(Mesh* parent_mesh)
: parent_mesh_(parent_mesh),
  local_part_(nullptr)
{
	this->init_element_vector(0);
	this->init_node_vector(0);
}


BCMesh::~BCMesh()
{
	if (local_part_!=nullptr) delete local_part_;
}


Range<ElementAccessor<3>> BCMesh::elements_range() const
{
	auto bgn_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(parent_mesh_, parent_mesh_->bulk_size_) );
	auto end_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(parent_mesh_, parent_mesh_->element_vec_.size()) );
    return Range<ElementAccessor<3>>(bgn_it, end_it);
}


unsigned int BCMesh::n_elements(FMT_UNUSED bool boundary) const {
	return parent_mesh_->element_ids_.size()-parent_mesh_->bulk_size_;
}


Partitioning *BCMesh::get_part() {
    return parent_mesh_->get_part();
}

const LongIdx *BCMesh::get_local_part() {
	if (local_part_ == nullptr) {
		local_part_ = new LongIdx[this->n_elements()];
		unsigned int bc_ele_idx;
		for (auto ele : parent_mesh_->elements_range())
			if (ele->boundary_idx_ != NULL)
				for (unsigned int i=0; i<ele->n_sides(); ++i)
					if ((int)ele->boundary_idx_[i] != -1) {
						bc_ele_idx = parent_mesh_->boundary_[ ele->boundary_idx_[i] ].bc_ele_idx_;
						local_part_[bc_ele_idx - parent_mesh_->bulk_size_] = parent_mesh_->get_local_part()[ele.idx()];
					}
	}
	return local_part_;
}


bool BCMesh::check_compatible_mesh( Mesh & mesh, vector<LongIdx> & bulk_elements_id, vector<LongIdx> & boundary_elements_id ) {
	return parent_mesh_->check_compatible_mesh(mesh, bulk_elements_id, boundary_elements_id);
}


unsigned int BCMesh::n_nodes() const {
    return parent_mesh_->nodes_.size();
}


ElementAccessor<3> BCMesh::element_accessor(unsigned int idx) const {
    return ElementAccessor<3>(parent_mesh_, idx+parent_mesh_->bulk_size_);
}


