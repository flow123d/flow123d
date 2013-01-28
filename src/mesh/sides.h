/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief ???
 *
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
    // Basic data
    //struct Boundary *cond; // Boundary condition  - if prescribed

    // Results
    //double flux; // Flux through side
    //double scalar; // Scalar quantity (piez. head or pressure)
    //double pscalar; // As scalar but in previous time step
    //struct Edge *edge_; // Edge to which belonged

    Side()
    : element_(NULL), el_idx_(0)
    {}

    inline Side(ElementIter ele, unsigned int set_lnum);
    double metric() const;
    arma::vec3 centre() const; // Centre of side
    arma::vec3 normal() const; // Vector of (generalized) normal

    inline unsigned int n_nodes() const;

    inline unsigned int dim() const;

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool is_external() const;

    inline const Node * node(unsigned int i) const;

    inline ElementFullIter element() const;

    inline Mesh * mesh() const;

    inline unsigned int edge_idx() const;

    inline Edge * edge() const;

    inline Boundary * cond() const;
    inline unsigned int cond_idx() const;

    inline unsigned int el_idx() const;

    inline bool valid() const;

    inline void inc();

    /// This is necessary by current DofHandler, should change this
    inline void *make_ptr() const;
private:

    arma::vec3 normal_point() const;
    arma::vec3 normal_line() const;
    arma::vec3 normal_triangle() const;

    // Topology of the mesh

    Element * element_; // Pointer to element to which belonged
    unsigned int el_idx_; // Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

    //Mesh    *mesh_;
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
