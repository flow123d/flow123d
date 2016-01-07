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

#include "mesh/nodes.hh"
#include "mesh/region.hh"
#include "mesh/bounding_box.hh"

template <int spacedim>
class ElementAccessor;

class Mesh;
class Side;
class SideIter;
class Neighbour;



//=============================================================================
// STRUCTURE OF THE ELEMENT OF THE MESH
//=============================================================================
class Element
{
public:
    Element();
    Element(unsigned int dim, Mesh *mesh_in, RegionIdx reg);
    void init(unsigned int dim, Mesh *mesh_in, RegionIdx reg);
    ~Element();


    inline unsigned int dim() const;
    inline unsigned int index() const;
    unsigned int n_sides() const; // Number of sides
    unsigned int n_nodes() const; // Number of nodes
    
    ///Gets ElementAccessor of this element
    ElementAccessor<3> element_accessor() const;
    
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
    double quality_measure_smooth();

    unsigned int n_sides_by_dim(unsigned int side_dim);
    inline SideIter side(const unsigned int loc_index);
    inline const SideIter side(const unsigned int loc_index) const;
    Region region() const;
    inline RegionIdx region_idx() const
        { return region_idx_; }
    
    unsigned int id() const;

    int pid; // Id # of mesh partition

    // Type specific data
    Node** node; // Element's nodes


    unsigned int *edge_idx_; // Edges on sides
    unsigned int *boundary_idx_; // Possible boundaries on sides (REMOVE) all bcd assembly should be done through iterating over boundaries
                           // ?? deal.ii has this not only boundary iterators
    /**
* Indices of permutations of nodes on sides.
* It determines, in which order to take the nodes of the side so as to obtain
* the same order as on the reference side (side 0 on the particular edge).
*
* Permutations are defined in RefElement::side_permutations.
*/
    unsigned int *permutation_idx_;

    /**
* Computes bounding box of element (OBSOLETE)
*/
    void get_bounding_box(BoundingBox &bounding_box) const;

    /**
* Return bounding box of the element.
*/
    inline BoundingBox bounding_box() {
     return BoundingBox(this->vertex_list());
    }

    /**
* Return list of element vertices.
*/
    inline vector<arma::vec3> vertex_list() {
     vector<arma::vec3> vertices(this->n_nodes());
     for(unsigned int i=0; i<n_nodes(); i++) vertices[i]=node[i]->point();
     return vertices;
    }


    unsigned int n_neighs_vb; // # of neighbours, V-B type (comp.)
                            // only ngh from this element to higher dimension edge
    Neighbour **neigh_vb; // List og neighbours, V-B type (comp.)


    Mesh *mesh_; // should be removed as soon as the element is also an Accessor


protected:
    // Data readed from mesh file
    RegionIdx region_idx_;
    unsigned int dim_;

    friend class GmshMeshReader;
    friend class Mesh;

    template<int spacedim, class Value>
    friend class Field;

};




#define FOR_ELEMENT_NODES(i,j) for((j)=0;(j)<(i)->n_nodes();(j)++)
#define FOR_ELEMENT_SIDES(i,j) for(unsigned int j=0; j < (i)->n_sides(); j++)
#define FOR_ELM_NEIGHS_VB(i,j) for((j)=0;(j)<(i)->n_neighs_vb;(j)++)


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
