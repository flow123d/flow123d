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
	this->nodes_ = parent_mesh_->nodes_;
	this->init_element_vector(1);
	this->init_node_vector(0);
}


BCMesh::~BCMesh()
{
	if (local_part_!=nullptr) delete local_part_;
}

void BCMesh::init_distribution()
{
	vector<LongIdx> loc_el_ids;

	// loop local parent elements and their sides to find boundary elements
	for (unsigned int iel = 0; iel<parent_mesh_->get_el_ds()->lsize(); iel++)
	{
		auto parent_el = parent_mesh_->element_accessor(parent_mesh_->get_el_4_loc()[iel]);
		for (unsigned int sid=0; sid<parent_el->n_sides(); sid++)
		{
			if (parent_el.side(sid)->is_boundary())
			{
				loc_el_ids.push_back(parent_el.side(sid)->cond().element_accessor().idx());
			}
		}
	}
	this->el_ds = new Distribution(loc_el_ids.size(), PETSC_COMM_WORLD);

	this->el_4_loc = new LongIdx[loc_el_ids.size()];
	for (unsigned int i=0; i<loc_el_ids.size(); i++) this->el_4_loc[i] = loc_el_ids[i];

	this->row_4_el = new LongIdx[n_elements()];
	vector<LongIdx> row_4_loc_el(n_elements(), 0);
	for (unsigned int i=0; i<loc_el_ids.size(); i++)
		row_4_loc_el[loc_el_ids[i]] = i + this->el_ds->begin();
	MPI_Allreduce(row_4_loc_el.data(), this->row_4_el, n_elements(), MPI_LONG_IDX, MPI_MAX, PETSC_COMM_WORLD);

}


Range<ElementAccessor<3>> BCMesh::elements_range() const
{
	auto bgn_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(parent_mesh_, 0, true) );
	auto end_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(parent_mesh_, element_vec_.size(), true) );
    return Range<ElementAccessor<3>>(bgn_it, end_it);
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
						local_part_[bc_ele_idx] = parent_mesh_->get_local_part()[ele.idx()];
					}
	}
	return local_part_;
}


std::shared_ptr<std::vector<LongIdx>> BCMesh::check_compatible_mesh( Mesh & input_mesh ) {
	return parent_mesh_->check_compatible_mesh(input_mesh);
}


std::shared_ptr<std::vector<LongIdx>> BCMesh::check_compatible_discont_mesh( Mesh & input_mesh) {
	return parent_mesh_->check_compatible_discont_mesh(input_mesh);
}


ElementAccessor<3> BCMesh::element_accessor(unsigned int idx) const {
    return ElementAccessor<3>(parent_mesh_, idx, true);
}


NodeAccessor<3> BCMesh::node(unsigned int) const
{
	ASSERT(	false );
	return NodeAccessor<3>();
}

Boundary BCMesh::boundary(unsigned int) const
{
	ASSERT( false );
	return Boundary();
}

const RegionDB &BCMesh::region_db() const
{
	ASSERT( false );
	static RegionDB r;
	return r;
}
