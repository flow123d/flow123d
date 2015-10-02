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
#include "mesh/mesh_types.hh"

//=============================================================================
// STRUCTURE OF THE SIDE OF THE MESH
//=============================================================================

class Mesh;

class Side {
public:
    Side()
    : element_(NULL), el_idx_(0)
    {}

    inline Side(const Element * ele, unsigned int set_lnum);
    double measure() const;
    arma::vec3 centre() const; // Centre of side
    arma::vec3 normal() const; // Vector of (generalized) normal

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
    inline const Node * node(unsigned int i) const;

    /**
     * Returns full iterator to the element of the side.
     */
    inline ElementFullIter element() const;

    /**
     * Returns pointer to the mesh.
     */
    inline Mesh * mesh() const;

    /**
     * Returns global index of the edge connected to the side.
     */
    inline unsigned int edge_idx() const;

    /**
     * Returns pointer to the edge connected to the side.
     */
    inline Edge * edge() const;

    inline Boundary * cond() const;
    inline unsigned int cond_idx() const;

    /**
     * Returns local index of the side on the element.
     */
    inline unsigned int el_idx() const;

    /**
     * Returns true if the side has assigned element.
     */
    inline bool valid() const;

    /**
     * Iterate over local sides of the element.
     */
    inline void inc();

    /// This is necessary by current DofHandler, should change this
    inline void *make_ptr() const;
private:

    arma::vec3 normal_point() const;
    arma::vec3 normal_line() const;
    arma::vec3 normal_triangle() const;

    // Topology of the mesh

    const Element * element_; // Pointer to element to which belonged
    unsigned int el_idx_; // Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

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

    inline bool operator==(const SideIter &other) {
        return (side_.element() == other.side_.element() ) && ( side_.el_idx() == other.side_.el_idx() );
    }


    inline bool operator!=(const SideIter &other) {
        return !( *this == other);
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
