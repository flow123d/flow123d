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

struct Side;
class Mesh;
struct Problem;

//=============================================================================
// STRUCTURE OF THE SIDE OF THE MESH
//=============================================================================



class Side {
public:
    // Basic data
    int id; // Id # of side
    int type; // INTERNAL | EXTERNAL
    int dim;
    // Topology of the mesh
    ElementIter element; // Pointer to element to which belonged
    int lnum; // Local # of side in element  (to remove it, we heve to remove calc_side_rhs)

    Node** node; // Pointers to sides's nodes

    struct Boundary *cond; // Boundary condition  - if prescribed
    struct Edge *edge; // Edge to wich belonged
    struct Neighbour *neigh_bv; // Neighbour, B-V type (comp.)
    // Geometry
    double normal[ 3 ]; // Vector of (generalized) normal
    // Matrix
    int c_row; // # of row in block C
    int c_col; // # of col in block C
    double c_val; // Value in block C
    // Results
    double flux; // Flux through side
    double scalar; // Scalar quantity (piez. head or pressure)
    double pscalar; // As scalar but in previous time step


    Side();
    double metric();
    arma::vec3 centre(); // Centre of side

    void reinit(ElementIter ele, unsigned int set_dim, int set_id, int set_lnum);
    inline unsigned int n_nodes() const
        {return dim+1;}


private:
    double length_line();
    double area_triangle();
};

#define EXTERNAL    0
#define INTERNAL    1



void make_side_list(Mesh*);
void side_calculation_mh(Mesh *mesh);
void calc_side_metrics(struct Side*);
void side_shape_specific(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
