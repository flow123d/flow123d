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
#include "mesh/intersection.hh"

#include "input/input_type.hh"
#include "input/accessors.hh"

// Forward declarations
template <int spacedim>
class ElementAccessor;
class GmshMeshReader;



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


class BoundarySegment {
public:
    static Input::Type::Record input_type;
};

//=============================================================================
// STRUCTURE OF THE MESH
//=============================================================================

class Mesh {
public:
    static Input::Type::Record input_type;

    Input::Record in_record_;

    /** Labels for coordinate indexes in arma::vec3 representing vectors and points.*/
    enum {x_coord=0, y_coord=1, z_coord=2};

    Mesh();
    Mesh(Input::Record in_record);
    void reinit(Input::Record in_record);

    inline unsigned int n_elements() const {
        return element.size();
    }

    inline unsigned int n_boundaries() const {
        return boundary.size();
    }

    inline unsigned int n_edges() const {
        return edge.size();
    }

    void read_intersections();
    void make_intersec_elements();
    // void make_edge_list_from_neigh();

    unsigned int n_sides();

    inline unsigned int n_vb_neighbours() const {
        return vb_neighbours_.size();
    }
    /**
     * Setup various links between mesh entities. Should be simplified.
     */
    void setup_topology(istream *in = NULL);

    /**
     *  This replaces read_neighbours() in order to avoid using NGH preprocessor.
     *
     *  TODO:
     *  - Avoid maps:
     *    1) node_elements can be vector, to get index from NodeIter call : mesh->node_vector.index(node_iter)
     *       in future we replace pointers by indexes, this is also safer for reallocations
     *    2) sort element vectors in node_elements
     *    3) make independent function that takes list of nodes and find intersection of their element lists:
     *       Mesh::intersect_element_lists(vector<unsigned int> &nodes_list, vector<unsigned int> &intersection_element_list)
     *       since element_lists are sorted we can find intersection easily using e.g std::set_intersection algorithm
     *       on the first pair of lists and then on remaining lists and resulting intersection list.
     *
     *    4) replace EdgeVector by std::vector<Edge>
     *
     *    5) need not to have temporary array for Edges, only postpone setting pointers in elements and set them
     *       after edges are found; we can temporary save Edge index instead of pointer in Neigbours and elements
     *
     *    6) Try replace Edge * by indexes in Neigbours and elements (anyway we have mesh pointer in elements so it is accessible also from Neigbours)
     *
     */
    void make_neighbours_and_edges();

    /**
     * This set pointers from elements to materials. Mesh should store only material IDs of indices.
     * This implies that element->volume can not be mesh property. Since fracture openning is material parameter.
     */
    void setup_materials( MaterialDatabase &base);
    void make_element_geometry();
    vector<int> const &all_elements_id();


    ElementAccessor<3> element_accessor(unsigned int idx, bool boundary=false);

    /// Vector of nodes of the mesh.
    NodeVector node_vector;
    /// Vector of elements of the mesh.
    ElementVector element;

    /// Vector of boundary sides where is prescribed boundary condition.
    /// TODO: apply all boundary conditions in the main assembling cycle over elements and remove this Vector.
    BoundaryVector boundary;
    /// vector of boundary elements - should replace 'boundary'
    /// TODO: put both bulk and bc elements (on zero level) to the same vector or make better map id->element for field inputs that use element IDs
    /// the avoid usage of ElementVector etc.
    ElementVector bc_elements;

    /// Vector of MH edges, this should not be part of the geometrical mesh
    EdgeVector edge;

    flow::VectorId<int> bcd_group_id; // gives a index of group for an id

    /**
     * Vector of individual intersections of two elements.
     * This is enough for local mortar.
     */
    vector<Intersection>  intersections;

    /**
     * For every element El we have vector of indices into @var intersections array for every intersection in which El is master element.
     * This is necessary for true mortar.
     */
    vector<vector<unsigned int> >  master_elements;

    vector<Neighbour> vb_neighbours_;
    int n_materials; // # of materials

    int n_insides; // # of internal sides
    int n_exsides; // # of external sides
    int n_sides_; // total number of sides (should be easy to count when we have separated dimensions

    int n_lines; // Number of line elements
    int n_triangles; // Number of triangle elements
    int n_tetrahedras; // Number of tetrahedra elements

    // for every side dimension D = 0 .. 2
    // for every element side 0 .. D+1
    // for every side node 0 .. D
    // index into element node array
    vector< vector< vector<unsigned int> > > side_nodes;
    
    string neigh_fname_;
    string bcd_fname_;

protected:

    void create_node_element_lists();
    void intersect_element_lists(vector<unsigned int> const &nodes_list, vector<unsigned int> &intersection_element_list);

    void element_to_neigh_vb();
    void create_external_boundary();

    void count_element_types();
    void count_side_types();

    unsigned int n_bb_neigh, n_vb_neigh;
    vector<Neighbour_both> neighbours_;

    /// Vector of both bulk and boundary IDs. Bulk elements come first, then boundary elements, but only the portion that appears
    /// in input mesh file and has ID assigned.
    ///
    /// TODO: Rather should be part of GMSH reader, but in such case we need store pointer to it in the mesh (good idea, but need more general interface for readers)
    vector<int> all_elements_id_;
    /// Number of elements read from input.
    unsigned int n_all_input_elements_;

    // For each node the vector contains a list of elements that use this node
    vector<vector<unsigned int> > node_elements;

    friend class GmshMeshReader;
};


#include "mesh/side_impl.hh"
#include "mesh/element_impls.hh"
#include "mesh/neighbours_impl.hh"

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


#define FOR_NEIGHBOURS(_mesh_, it) \
    for( std::vector<Neighbour>::iterator it = _mesh_->vb_neighbours_.begin(); \
         (it)!= _mesh_->vb_neighbours_.end(); ++it)

#define FOR_NEIGH_ELEMENTS(i,j) for((j)=0;(j)<(i)->n_elements;(j)++)
#define FOR_NEIGH_SIDES(i,j)    for((j)=0;(j)<(i)->n_sides;(j)++)

//int *max_entry();


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
