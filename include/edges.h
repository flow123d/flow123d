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

#ifndef MAKE_EDGES_H
#define MAKE_EDGES_H

#include "mesh.h"

//=============================================================================
// STRUCTURE OF THE EDGE OF THE MESH
//=============================================================================
typedef struct Edge
{
    // Basic
    int  id;        // Id # of the edge
    // Topology of the mesh
    int  n_sides;   // # of sides of edge
    struct Side **side; // sides of edge (could be more then two e.g. 1D mesh in 2d space with crossing )
    struct Neighbour *neigh_vb; // "Compatible" neighbouring
    struct Neighbour *neigh_bb; // ??? this is what
    // List
    struct Edge *prev;  // Previous edge in the list
    struct Edge *next;  // Next edge in the list
    // Matrix
    int  c_row;     // # of row in block C (and E and F) (MH)
    double  f_val;      // diagonal value  in block F
    double  f_rhs;      // rhs value
    // Misc
    int      aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
} Edge;

#define FOR_EDGES(i)        for((i)=mesh->edge;(i)!=NULL;(i)=(i)->next)
#define FOR_EDGE_SIDES(i,j) for((j)=0;(j)<(i)->n_sides;(j)++)

void make_edge_list(Mesh*);
void edge_calculation_mh(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
