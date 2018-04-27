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
 * @file    node_accessor.hh
 * @brief
 */

#ifndef NODE_ACCESSOR_HH_
#define NODE_ACCESSOR_HH_

#include "mesh/nodes.hh"
#include "mesh/mesh.h"

/**
 * Node accessor templated just by dimension of the embedding space, used for access to nodes out of Mesh.
 * This should allow algorithms over nodes.
 */
template <int spacedim>
class NodeAccessor {
public:
    /**
     * Default invalid accessor.
     */
	NodeAccessor()
    : mesh_(NULL)
    {}

    /**
     * Node accessor.
     */
	NodeAccessor(const Mesh *mesh, unsigned int idx)
    : mesh_(mesh), node_idx_(idx)
    {}

    inline bool is_valid() const {
        return mesh_ != NULL;
    }

    inline const Node * node() const {
        return &(mesh_->node_vec_[node_idx_]);
    }

    inline unsigned int idx() const {
        return node_idx_;
    }

    inline unsigned int index() const {
    	return (unsigned int)mesh_->find_node_id(node_idx_);
    }

    inline void inc() {
        ASSERT(!is_valid()).error("Do not call inc() for invalid accessor!");
        node_idx_++;
    }

    bool operator==(const NodeAccessor<spacedim>& other) {
    	return (node_idx_ == other.node_idx_);
    }

    /**
     * -> dereference operator
     *
     * Allow simplify calling of node() method. Example:
 @code
     NodeAccessor<3> node_ac(mesh, index);
     arma::vec centre;
     centre = node_ac.node()->point();  // full format of access to node
     centre = node_ac->point();         // short format with dereference operator
 @endcode
     */
    inline const Node * operator ->() const {
    	return &(mesh_->node_vec_[node_idx_]);
    }

private:
    /// Pointer to the mesh owning the node.
    const Mesh *mesh_;

    /// Index into Mesh::node_vec_ array.
    unsigned int node_idx_;
};


#endif /* NODE_ACCESSOR_HH_ */
