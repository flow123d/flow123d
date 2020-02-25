/*!
 *
 * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    edge.hh
 * @brief   
 */

#ifndef MAKE_EDGES_H
#define MAKE_EDGES_H

#include "mesh/sides.h"     // for SideIter
#include "mesh/mesh.h"      // for EdgeData

//=============================================================================
// STRUCTURE OF THE EDGE OF THE MESH
//=============================================================================
class Edge
{
public:
    Edge()
    : mesh_(nullptr),
      edge_idx_(Mesh::undef_idx)
    {}

    Edge(const Mesh *mesh, unsigned int edge_idx)
    : mesh_(mesh),
      edge_idx_(edge_idx)
    {}

    inline bool is_valid() const {
        return mesh_ != nullptr;
    }

    inline unsigned int idx() const {
        ASSERT_DBG(is_valid());
        return edge_idx_;
    }

    inline unsigned int n_sides() const
    { return edge_data()->n_sides;}

    inline SideIter side(const unsigned int i) const {
        return edge_data()->side_[i];
    }

    inline void inc() {
        ASSERT(is_valid()).error("Do not call inc() for invalid accessor!");
        edge_idx_++;
    }

    bool operator==(const Edge& other) const{
    	return (edge_idx_ == other.edge_idx_);
    }

private:
    /// Pointer to the mesh owning the node.
    const Mesh *mesh_;
    /// Index into Mesh::edges vector.
    unsigned int edge_idx_;

    inline const EdgeData* edge_data() const
    {
        ASSERT_DBG(is_valid());
        ASSERT_LT_DBG(edge_idx_, mesh_->edges.size());
        return &mesh_->edges[edge_idx_];
    }
};


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
