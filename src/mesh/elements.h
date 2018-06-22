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
 * @file    elements.h
 * @brief   
*/

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <ext/alloc_traits.h>                  // for __alloc_traits<>::valu...
#include <string.h>                            // for memcpy
#include <new>                                 // for operator new[]
#include <ostream>                             // for operator<<
#include <string>                              // for operator<<
#include <vector>                              // for vector
#include <armadillo>
#include "side_impl.hh"
#include "mesh/bounding_box.hh"                // for BoundingBox
#include "mesh/nodes.hh"                       // for Node
#include "mesh/ref_element.hh"                 // for RefElement
#include "mesh/region.hh"                      // for RegionIdx, Region
#include "system/asserts.hh"                   // for Assert, ASSERT
#include "sides.h"
#include "mesh/mesh.h"

//class Mesh;
class Neighbour;
//class SideIter;
template <int spacedim> class ElementAccessor;



//=============================================================================
// STRUCTURE OF THE ELEMENT OF THE MESH
//=============================================================================
class Element
{
public:
    Element();
    Element(unsigned int dim, RegionIdx reg);
    void init(unsigned int dim, RegionIdx reg);
    ~Element();


    inline unsigned int dim() const;
    unsigned int n_sides() const; // Number of sides
    unsigned int n_nodes() const; // Number of nodes
    
    /// Computes the measure of the element.
    double measure() const;
    
    /** Computes the Jacobian of the element.
     * J = det ( 1  1  1  1 )
     *           x1 x2 x3 x4
     *           y1 y2 y3 y4
     *           z1 z2 z3 z4
     */
    double tetrahedron_jacobian() const;
    
    /// Computes the barycenter.
    arma::vec3 centre() const;
    /**
* Quality of the element based on the smooth and scale-invariant quality measures proposed in:
* J. R. Schewchuk: What is a Good Linear Element?
*
* We scale the measure so that is gives value 1 for regular elements. Line 1d elements
* have always quality 1.
*/
    double quality_measure_smooth(SideIter side) const;

    inline RegionIdx region_idx() const
        { return region_idx_; }
    
    /**
     * Computes bounding box of element (OBSOLETE) ??
     */
    void get_bounding_box(BoundingBox &bounding_box) const;

    /// Return precomputed bounding box.
    //BoundingBox &get_bounding_box_fast(BoundingBox &bounding_box) const;

    /**
    * Return bounding box of the element.
    * Simpler code, but need to check performance penelty.
    */
    inline BoundingBox bounding_box() const {
     return BoundingBox(this->vertex_list());
    }

    /**
     * Return list of element vertices.
     */
    inline vector<arma::vec3> vertex_list() const {
    	vector<arma::vec3> vertices(this->n_nodes());
    	for(unsigned int i=0; i<n_nodes(); i++) vertices[i]=node[i]->point();
    	return vertices;
    }
    
    /// Return edge_idx of given index
    inline unsigned int edge_idx(unsigned int edg_idx) const;

    /// Return permutation_idx of given index
    inline unsigned int permutation_idx(unsigned int prm_idx) const;

    /// Return Id of mesh partition
    inline int pid() const {
    	return pid_;
    }

    /// Return number of neighbours
    inline unsigned int n_neighs_vb() const {
    	return n_neighs_vb_;
    }

    //unsigned int get_proc() const;


    // TODO move data members to protected part, add access trough getters or use direct access of friend class Mesh

    // Type specific data
    Node** node; // Element's nodes
        // TODO remove direct access in many places in fem, fields, flow, intersection, io, mesh, transport and unit tests


    unsigned int *boundary_idx_; // Possible boundaries on sides (REMOVE) all bcd assembly should be done through iterating over boundaries
                           // ?? deal.ii has this not only boundary iterators
                           // TODO remove direct access in balance, side and transport

    Neighbour **neigh_vb; // List og neighbours, V-B type (comp.)
        // TODO remove direct access in DarcyFlow, MhDofHandler, Mesh? Partitioning and Trabsport


protected:
    int pid_;                            ///< Id # of mesh partition
    std::vector<unsigned int> edge_idx_; ///< Edges on sides
    mutable unsigned int n_neighs_vb_;   ///< # of neighbours, V-B type (comp.)
                                         // only ngh from this element to higher dimension edge
                                         // TODO fix and remove mutable directive

    /**
    * Indices of permutations of nodes on sides.
    * It determines, in which order to take the nodes of the side so as to obtain
    * the same order as on the reference side (side 0 on the particular edge).
    *
    * Permutations are defined in RefElement::side_permutations.
    *
    * TODO fix and remove mutable directive
    */
    mutable std::vector<unsigned int> permutation_idx_;

    // Data readed from mesh file
    RegionIdx  region_idx_;
    unsigned int dim_;

    friend class Mesh;

    template<int spacedim, class Value>
    friend class Field;

};


inline unsigned int Element::dim() const {
    return dim_;
}


inline unsigned int Element::n_nodes() const {
    return dim()+1;
}



inline unsigned int Element::n_sides() const {
    return dim()+1;
}

inline unsigned int Element::edge_idx(unsigned int edg_idx) const {
	ASSERT(edg_idx<edge_idx_.size())(edg_idx)(edge_idx_.size()).error("Index of Edge is out of bound!");
	return edge_idx_[edg_idx];
}

inline unsigned int Element::permutation_idx(unsigned int prm_idx) const {
	ASSERT(prm_idx<permutation_idx_.size())(prm_idx)(permutation_idx_.size()).error("Index of permutation is out of bound!");
	return permutation_idx_[prm_idx];
}

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
