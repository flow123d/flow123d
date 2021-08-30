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

	make_neighbours_and_edges();
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


void BCMesh::make_neighbours_and_edges()
{
    Neighbour neighbour;
    EdgeData *edg = nullptr;
    unsigned int ngh_element_idx;
    unsigned int last_edge_idx = Mesh::undef_idx;

    neighbour.mesh_ = this;

    create_node_element_lists();

	// pointers to created edges
    edges.resize(0); // be sure that edges are empty

	vector<unsigned int> side_nodes;
	vector<unsigned int> intersection_list; // list of elements in intersection of node element lists

	// Now we go through all element sides and create edges and neighbours
	for (auto e : this->elements_range()) {
		for (unsigned int s=0; s<e->n_sides(); s++)
		{
			// skip sides that were already found
			if (e->edge_idx(s) != Mesh::undef_idx) continue;

			// Find all elements that share this side.
			side_nodes.resize(e.side(s)->n_nodes());
			for (unsigned n=0; n<e.side(s)->n_nodes(); n++) side_nodes[n] = e.side(s)->node(n).idx();
			intersect_element_lists(side_nodes, intersection_list);

			bool is_neighbour = find_lower_dim_element(intersection_list, e->dim(), ngh_element_idx);

			if (is_neighbour) { // edge connects elements of different dimensions
				// Initialize for the neighbour case.
			    neighbour.elem_idx_ = ngh_element_idx;
            } else { // edge connects only elements of the same dimension
                // Initialize for the edge case.
                last_edge_idx=edges.size();
                edges.resize(last_edge_idx+1);
                edg = &( edges.back() );
                edg->n_sides = 0;
                edg->side_ = new struct SideIter[ intersection_list.size() ];
                if (e->dim() > 0)
					if (intersection_list.size() > max_edge_sides_[e->dim()-1])
                		max_edge_sides_[e->dim()-1] = intersection_list.size();

                if (intersection_list.size() <= 1) {
                	// outer edge, create boundary object as well
                    edg->n_sides=1;
                    edg->side_[0] = e.side(s);
                    element_vec_[e.idx()].edge_idx_[s] = last_edge_idx;

                    continue; // next side of element e
                }
			}

			// go through the elements connected to the edge or neighbour
			// setup neigbour or edge
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect) {
            	ElementAccessor<3> elem = this->element_accessor(*isect);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    if (elem->edge_idx(ecs) != Mesh::undef_idx) continue; // ??? This should not happen.
                    SideIter si = elem.side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (is_neighbour) {
                            // create a new edge and neighbour for this side, and element to the edge
                            last_edge_idx=edges.size();
                            edges.resize(last_edge_idx+1);
                            edg = &( edges.back() );
                            edg->n_sides = 1;
                            edg->side_ = new struct SideIter[1];
                            edg->side_[0] = si;
                            element_vec_[elem.idx()].edge_idx_[ecs] = last_edge_idx;

                            neighbour.edge_idx_ = last_edge_idx;

                            vb_neighbours_.push_back(neighbour); // copy neighbour with this edge setting
                        } else {
                            // connect the side to the edge, and side to the edge
                            ASSERT_PTR_DBG(edg);
                            edg->side_[ edg->n_sides++ ] = si;
                            ASSERT_DBG(last_edge_idx != Mesh::undef_idx);
                            element_vec_[elem.idx()].edge_idx_[ecs] = last_edge_idx;
                        }
                        break; // next element from intersection list
                    }
                } // search for side of other connected element
            } // connected elements

            if (! is_neighbour)
				ASSERT_EQ( (unsigned int) edg->n_sides, intersection_list.size())(e.index())(s).error("Missing edge sides.");
		} // for element sides
	}   // for elements

	MessageOut().fmt( "Created {} edges and {} neighbours on boundary mesh.\n", edges.size(), vb_neighbours_.size() );
}
