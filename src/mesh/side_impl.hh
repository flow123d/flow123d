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
 * @file    side_impl.hh
 * @brief   
 */

#ifndef SIDE_IMPL_HH_
#define SIDE_IMPL_HH_

#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/boundaries.h"


inline Side::Side(const Mesh * mesh, unsigned int elem_idx, unsigned int set_lnum)
: mesh_(mesh), elem_idx_(elem_idx), side_idx_(set_lnum)
{
	mesh_->check_element_size(elem_idx);
}


    inline unsigned int Side::n_nodes() const {
        return dim()+1;
    }

    inline unsigned int Side::dim() const {
        return element()->dim()-1;
    }

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool Side::is_external() const {
        return edge().n_sides() == 1;
    }

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool Side::is_boundary() const {
        return is_external() && cond_idx() != Mesh::undef_idx;
    }

    inline NodeAccessor<3> Side::node(unsigned int i) const {
        int i_n = mesh_->side_nodes[dim()][side_idx_][i];

        return element().node_accessor( i_n );
    }

    inline ElementAccessor<3> Side::element() const {
    	ASSERT( valid() ).error("Wrong use of uninitialized accessor.\n");
        return mesh_->element_accessor( elem_idx_ );
    }

    inline const Mesh * Side::mesh() const {
        return this->mesh_;
    }

    inline unsigned int Side::edge_idx() const {
        return element()->edge_idx(side_idx_);
    }

    // inline const Edge * Side::edge() const {
    //     if (edge_idx() == Mesh::undef_idx) return NULL;
    //     else return &( mesh_->edges[ edge_idx() ] );
    // }

    // inline Boundary * Side::cond() const {
    //         if (cond_idx() == Mesh::undef_idx) return NULL;
    //         else return &( mesh_->boundary_[ cond_idx() ] );
    // }

    inline Edge Side::edge() const {
        return mesh_->edge(edge_idx());
    }

    inline Boundary Side::cond() const {
        return mesh_->boundary(cond_idx());
    }

    inline unsigned int Side::cond_idx() const {
         if (element()->boundary_idx_ == NULL) return Mesh::undef_idx;
         else return element()->boundary_idx_[side_idx_];
    }


    inline unsigned int Side::side_idx() const {
        return side_idx_;
    }


    inline unsigned int Side::elem_idx() const {
        return elem_idx_;
    }


    inline bool Side::valid() const {
        return mesh_!= NULL;
    }


    inline void Side::inc() {
    	side_idx_++;
    }


    //inline void *Side::make_ptr() const {
    //    return (void *)((long int) element() << (2 + side_idx_) );
    //}
#endif /* SIDE_IMPL_HH_ */
