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

class Mesh;
class Edge;

#include <mesh_types.hh>

//=============================================================================
// STRUCTURE OF THE SIDE OF THE MESH
//=============================================================================



class Side {
public:
    // Basic data
    //int id; // Id # of side
    //int type; // INTERNAL | EXTERNAL


    //Node** node; // Pointers to sides's nodes

    struct Boundary *cond; // Boundary condition  - if prescribed

    // Results
    //double flux; // Flux through side
    //double scalar; // Scalar quantity (piez. head or pressure)
    //double pscalar; // As scalar but in previous time step
    struct Edge *edge_; // Edge to which belonged


    Side();
    double metric() const;
    arma::vec3 centre() const; // Centre of side
    arma::vec3 normal() const; // Vector of (generalized) normal

    void reinit(Mesh *mesh, ElementIter ele,  int set_id, int set_lnum);

    inline unsigned int n_nodes() const;

    inline unsigned int dim() const;

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool is_external() const;

    inline const Node * node(unsigned int i) const;

    inline ElementFullIter element(); // unfortunately we can not have const here, since there are plenty of ELEMENT_FULL_ITER

    inline const Mesh * mesh() const;

    inline Edge * edge();  // unfortunately we can not have const here, we need to convert it to Full iterator

    inline unsigned int el_idx() const;


private:

    arma::vec3 normal_point() const;
    arma::vec3 normal_line() const;
    arma::vec3 normal_triangle() const;

    // Topology of the mesh

    Element * element_; // Pointer to element to which belonged
    unsigned int el_idx_; // Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

    Mesh    *mesh_;
};

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
