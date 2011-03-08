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

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "nodes.h"
#include <materials.hh>

struct Problem;
class Mesh;
struct MaterialDatabase;

//=============================================================================
// STRUCTURE OF THE ELEMENT OF THE MESH
//=============================================================================
typedef struct Element
{
    // Data readed from mesh file
    int      type;      //
    int      mid;       // Id # of material
    int      rid;       // Id # of region
    int      pid;       // Id # of mesh partition
   // double   size;      // Area of 1D elm or height of 2D elm

    // Type specific data
    int dim;        // 1 or 2 or 3
    int n_sides;    // Number of sides

    int    n_nodes; // Number of nodes
    Node** node;    // Element's nodes

    MaterialDatabase::Iter material; // Element's material
    struct Side **side; // Element's sides
    int      n_neighs_vv;   // # of neighbours, V-V type (noncomp.)
    struct Neighbour **neigh_vv; // List og neighbours, V-V type (noncomp.)
    int      n_neighs_vb;   // # of neighbours, V-B type (comp.)
                            // only ngh from this element to higher dimension edge
    struct Neighbour **neigh_vb; // List og neighbours, V-B type (comp.)
    struct Concentration  *start_conc;    // Start concentration  - if prescribed
    int n_subst;        // Number of substances
    // Material properties
    double   *k;    // Tensor of hydraulic conductivity
    double   *a;    // Tensor of hydraulic resistance (k^(-1))
                    // --------------------------
                    // matrices a,k,loc are stored as vectors, and should be retyped
                    // if one needs double index access, see SmallMtxN types in math_fce.h
    double   stor;          // Storativity
    // Geometrical properties
    double   measure;   // Metrics of the element (length, area, volume)
    double   volume;        // Volume of the element
    double   centre[ 3 ];   // Centre of mass of element
    // Parameters of the basis functions
    double   *bas_alfa;      // Parameters alfa
    double   *bas_beta;      // Parameters beta
    double   *bas_gama;      // Parameters gama
    double   *bas_delta;      // Parameters delta
    // Matrix
    double *loc;        // Local matrix
    double *loc_inv;    // Inverse of the local matrix
    double  *rhs;       // Rhs - vector q1, gradients
    //double   rhs_b;     // Rhs - vector q2, sources
    // int      rhs_b_stat_index;   // Index of the item in the matrix where value is stored (Matrix->rhs)
    // double   rhs_b_stat;    // Value of rhs_b before unsteady flow changes
    int  a_row;     // # of first row in the block A
    int      b_row;     // # of row in the block B
    // int      b_col;     // # of first column in the block B
    // double  *b_val;     // Values in block B
    int  d_row_count;   // # of nonzeros in row of the block D
    int *d_col;     // columns with nonzeros in D
    double  *d_val; // values of entries in D
    int *d_el;      // ids of noncompatible conected elements (itself reference included)

    // int      d_val_stat_index;  // Index of the item in the matrix where value is stored (Matrix->a)
    // double   d_val_stat;        // value of diag. entry in D; computed as steady flow value
    // int      d_diag;            // Index of value in the array d_val where entry on the diagonal is stored
    int  e_row_count;   // # of nonzeros in row of the block E
    int *e_col;     // columns with nonzeros in E
    int *e_edge_id;   // ids of compatible conected edeges
    double  *e_val;     // values of entries in E
    // Time terms - diagonal in D block and additional RHS
    double tAddRHS, tAddDiag;
    
    // Results
    double   vector[ 3 ];   // Vector quantity - velocity of flow
    double   v_length;  // Length of vector quantity
    double   scalar;    // Scalar quantity (piez. head or pressure)
    double   pscalar;   // As scalar but in previous time step
    double   balance;   // Balance of flux

    double   scalar_it; // for density iteration
    double   *conc_prev;         // for density iteration
    double   *conc_prev_immobile;         // for density iteration
    double   *conc_prev_sorb;         // for density iteration
    double   *conc_prev_immobile_sorb;         // for density iteration

    // Misc
    int     aux;       // Auxiliary flag
    double  faux;      // Auxiliary number
} Element;


#define D_DIAG 0        // in D block diagonal is always at zero position in d_val,d_col

#define xPOINT 0
#define LINE 1
#define TRIANGLE 2
#define TETRAHEDRON 4

#define PROP_S 1        //Area of 1D element
#define PROP_H 2        //Height of 2D element
#define PROP_V 3        //Volume of 3D element

#define FOR_ELEMENT_NODES(i,j)  for((j)=0;(j)<(i)->n_nodes;(j)++)
#define FOR_ELEMENT_SIDES(i,j)  for(unsigned int j=0; j < (i)->n_sides; j++)
#define FOR_ELM_NEIGHS_VV(i,j)  for((j)=0;(j)<(i)->n_neighs_vv;(j)++)
#define FOR_ELM_NEIGHS_VB(i,j)  for((j)=0;(j)<(i)->n_neighs_vb;(j)++)

void read_element_list(Mesh*);
void element_calculation_mh(Mesh*);
void make_element_geometry();
void element_calculation_unsteady(struct Problem*);

//void read_element_properties(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
