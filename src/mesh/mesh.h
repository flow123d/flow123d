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

#ifndef MAKE_MESH_H
#define MAKE_MESH_H

#include <vector>
#include "mesh/mesh_types.hh"

#include "mesh/nodes.hh"
#include "mesh/elements.h"
#include "mesh/sides.h"
#include "mesh/edges.h"
#include "mesh/neighbours.h"
#include "mesh/boundaries.h"


#define ELM  0
#define BC  1
#define NODE  2

/**
 *  This parameter limits volume of elements from below.
 */
#define MESH_CRITICAL_VOLUME 1.0E-12

/**
 * Provides for statement to iterate over the Nodes of the Mesh.
 * The parameter is FullIter local variable of the cycle, so it need not be declared before.
 * Macro assume that variable Mesh *mesh; is declared and points to a valid Mesh structure.
 */
#define FOR_NODES(_mesh_, i) \
    for( NodeFullIter i( _mesh_->node_vector.begin() ); \
        i != _mesh_->node_vector.end(); \
        ++i)

/**
 * Macro for conversion form Iter to FullIter for nodes.
 */
#define NODE_FULL_ITER(_mesh_,i) \
    _mesh_->node_vector.full_iter(i)

/**
 * Macro to get "NULL" ElementFullIter.
 */
#define NODE_FULL_ITER_NULL(_mesh_) \
    NodeFullIter(_mesh_->node_vector)

/**
 * Macro for conversion form Iter to FullIter for elements.
 */
#define ELEM_FULL_ITER(_mesh_,i) \
    _mesh_->element.full_iter(i)


#define FOR_NODE_ELEMENTS(i,j)   for((j)=0;(j)<(i)->n_elements();(j)++)
#define FOR_NODE_SIDES(i,j)      for((j)=0;(j)<(i)->n_sides;(j)++)

//=============================================================================
// STRUCTURE OF THE MESH
//=============================================================================

class Mesh {
private:

public:
    /** Labels for coordinate indexes in arma::vec3 representing vectors and points.*/
    enum {x_coord=0, y_coord=1, z_coord=2};

    Mesh();

    inline unsigned int n_elements() const {
        return element.size();
    }

    inline unsigned int n_boundaries() const {
        return boundary.size();
    }

    inline unsigned int n_edges() const {
        return edge.size();
    }

    unsigned int n_sides();
    /**
     * Setup various links between mesh entities. Should be simplified.
     */
    void setup_topology();

    /**
     * This set pointers from elements to materials. Mesh should store only material IDs of indices.
     * This implies that element->volume can not be mesh property. Since fracture openning is material parameter.
     */
    void setup_materials( MaterialDatabase &base);
    void make_element_geometry();

    // Files
    // DF - Move to ConstantDB
    // char *geometry_fname; // Name of file of nodes and elems
    // char *concentration_fname;//Name of file of concentration
    // char *transport_bcd_fname;//Name of file of transport BCD

    /// Vector of nodes of the mesh.
    NodeVector node_vector;
    /// Vector of elements of the mesh.
    ElementVector element;
    /// Vector of boundary sides where is prescribed boundary condition.
    /// TODO: apply all boundary conditions in the main assembling cycle over elements and remove this Vector.
    BoundaryVector boundary;
    /// Vector of MH edges, this should not be part of the geometrical mesh
    EdgeVector edge;

    flow::VectorId<int> bcd_group_id; // gives a index of group for an id

    int n_materials; // # of materials

    int n_insides; // # of internal sides
    int n_exsides; // # of external sides
    int n_sides_; // total number of sides (should be easy to count when we have separated dimensions

    int n_neighs;
    struct Neighbour *neighbour; // First neighbour
    struct Neighbour *l_neighbour; // Last neighbour

    int n_lines; // Number of line elements
    int n_triangles; // Number of triangle elements
    int n_tetrahedras; // Number of tetrahedra elements

    // for every side dimension D = 0 .. 2
    // for every element side 0 .. D+1
    // for every side node 0 .. D
    // index into element node array
    vector< vector< vector<unsigned int> > > side_nodes;

private:
    void count_element_types();

};


#include "mesh/side_impl.hh"
#include "element_impls.hh"

/**
 * Provides for statement to iterate over the Elements of the Mesh.
 * The parameter is FullIter local variable of the cycle, so it need not be declared before.
 * Macro assume that variable Mesh *mesh; is declared and points to a valid Mesh structure.
 */
#define FOR_ELEMENTS(_mesh_,__i) \
    for( ElementFullIter __i( _mesh_->element.begin() ); \
        __i != _mesh_->element.end(); \
        ++__i)

/**
 * Macro for conversion form Iter to FullIter for elements.
 */
#define ELEMENT_FULL_ITER(_mesh_,i) \
    _mesh_->element.full_iter(i)

/**
 * Macro to get "NULL" ElementFullIter.
 */
#define ELEMENT_FULL_ITER_NULL(_mesh_) \
    ElementFullIter(_mesh_->element)


#define FOR_BOUNDARIES(_mesh_,i) \
for( BoundaryFullIter i( _mesh_->boundary.begin() ); \
    i != _mesh_->boundary.end(); \
    ++i)

/**
 * Macro for conversion form Iter to FullIter for boundaries.
 */
#define BOUNDARY_FULL_ITER(_mesh_,i) \
    _mesh_->boundary.full_iter(i)

/**
 * Macro to get "NULL" BoundaryFullIter.
 */
#define BOUNDARY_NULL(_mesh_) \
    BoundaryFullIter(_mesh_->boundary)


/**
 * Provides for statement to iterate over the Edges of the Mesh. see FOR_ELEMENTS
 */
#define FOR_EDGES(_mesh_,__i) \
    for( EdgeFullIter __i( _mesh_->edge.begin() ); \
        __i !=_mesh_->edge.end(); \
        ++__i)

#define FOR_SIDES(_mesh_, it) \
    FOR_ELEMENTS(_mesh_, ele)  \
        for(SideIter it = ele->side(0); it->el_idx() < ele->n_sides(); ++it)

#define FOR_SIDE_NODES(i,j) for((j)=0;(j)<(i)->n_nodes;(j)++)


#define FOR_NEIGHBOURS(_mesh_, i)   for((i)=_mesh_->neighbour;(i)!=NULL;(i)=(i)->next)
#define FOR_NEIGH_ELEMENTS(i,j) for((j)=0;(j)<(i)->n_elements;(j)++)
#define FOR_NEIGH_SIDES(i,j)    for((j)=0;(j)<(i)->n_sides;(j)++)

//int *max_entry();


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
