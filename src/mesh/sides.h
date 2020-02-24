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
 * @file    sides.h
 * @brief   
 */

#ifndef SIDES_H
#define SIDES_H


#include <armadillo>
#include <stddef.h>                          // for NULL
//#include "mesh/mesh.h"
//#include "mesh/accessors.hh"

class Boundary;
class Edge;
class Element;
template <int spacedim> class NodeAccessor;
template <int spacedim> class ElementAccessor;

//=============================================================================
// STRUCTURE OF THE SIDE OF THE MESH
//=============================================================================

class Mesh;

class Side {
public:
    Side()
    : mesh_(NULL), elem_idx_(0), side_idx_(0)
    {}

    inline Side(const Mesh * mesh, unsigned int elem_idx, unsigned int set_lnum); ///< Constructor
    double measure() const;    ///< Calculate metrics of the side
    arma::vec3 centre() const; ///< Centre of side
    arma::vec3 normal() const; ///< Vector of (generalized) normal
    double diameter() const;   ///< Calculate the side diameter.

    /**
     * Returns number of nodes of the side.
     */
    inline unsigned int n_nodes() const;

    /**
     * Returns dimension of the side, that is dimension of the element minus one.
     */
    inline unsigned int dim() const;

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool is_external() const;

    /**
     * Returns node for given local index @p i on the side.
     */
    inline NodeAccessor<3> node(unsigned int i) const;

    /**
     * Returns iterator to the element of the side.
     */
    inline ElementAccessor<3> element() const;

    /**
     * Returns pointer to the mesh.
     */
    inline const Mesh * mesh() const;

    /**
     * Returns global index of the edge connected to the side.
     */
    inline unsigned int edge_idx() const;

    /**
     * Returns pointer to the edge connected to the side.
     */
    inline const Edge * edge() const;

    inline Boundary * cond() const;
    inline unsigned int cond_idx() const;

    /**
     * Returns local index of the side on the element.
     */
    inline unsigned int side_idx() const;

    /**
     * Returns index of element in Mesh::element_vec_.
     */
    inline unsigned int elem_idx() const;

    /**
     * Returns true if the side has assigned element.
     */
    inline bool valid() const;

    /**
     * Iterate over local sides of the element.
     */
    inline void inc();

    inline bool operator==(const Side &other) const {
        return (mesh_ == other.mesh_ ) && ( elem_idx_ == other.elem_idx_ )
        		&& ( side_idx_ == other.side_idx_ );
    }

    inline bool operator!=(const Side &other) const {
        return !( *this == other);
    }

    /// This is necessary by current DofHandler, should change this
    //inline void *make_ptr() const;
private:

    arma::vec3 normal_point() const;
    arma::vec3 normal_line() const;
    arma::vec3 normal_triangle() const;

    // Topology of the mesh

    const Mesh * mesh_;     ///< Pointer to Mesh to which belonged
    unsigned int elem_idx_; ///< Index of element in Mesh::element_vec_
    unsigned int side_idx_; ///< Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

};


/*
 * Iterator to a side.
 */
class SideIter {
public:
    SideIter()
    {}

    inline SideIter(const Side &side)
    : side_(side)
    {}

    inline bool operator==(const SideIter &other) const {
        return side_ == other.side_;
    }


    inline bool operator!=(const SideIter &other) const {
        return side_ != other.side_;
    }

    ///  * dereference operator
    inline const Side & operator *() const
            { return side_; }

    /// -> dereference operator
    inline const Side * operator ->() const
            { return &side_; }

    /// prefix increment iterate only on local element
    inline SideIter &operator ++ () {
        side_.inc();
        return (*this);
    }

private:
    Side side_;
};
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
