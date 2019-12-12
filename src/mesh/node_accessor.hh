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

#include "system/armor.hh"
#include "mesh/point.hh"
#include "mesh/mesh.h"


/**
 * Node accessor templated just by dimension of the embedding space, used for access to nodes out of Mesh.
 * This should allow algorithms over nodes.
 */
template <int spacedim>
class NodeAccessor {
public:
    typedef typename Space<spacedim>::Point Point;

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

    inline unsigned int idx() const {
        return node_idx_;
    }

    // TODO: rename to gmsh_id
    inline unsigned int index() const {
    	return (unsigned int)mesh_->find_node_id(node_idx_);
    }

    inline void inc() {
        ASSERT(is_valid()).error("Do not call inc() for invalid accessor!");
        node_idx_++;
    }

    bool operator==(const NodeAccessor<spacedim>& other) {
    	return (node_idx_ == other.node_idx_);
    }

    inline Point operator*() const
    { return mesh_->nodes_.vec<spacedim>(node_idx_); }


private:
    /// Pointer to the mesh owning the node.
    const Mesh *mesh_;

    /// Index into Mesh::node_vec_ array.
    unsigned int node_idx_;
};


#endif /* NODE_ACCESSOR_HH_ */
