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
 * @brief    Various element oriented stuff, should be restricted to purely geometric functions
 * @ingroup mesh
 *
 */

#include <vector>

#include <string>

#include "system/system.hh"
#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "elements.h"
#include "element_impls.hh"

// following deps. should be removed
#include "mesh/boundaries.h"
#include "materials.hh"
#include "mesh/accessors.hh"



Element::Element()
:  pid(0),

  node(NULL),

  material(NULL),
  edge_idx_(NULL),
  boundary_idx_(NULL),

  n_neighs_vb(0),
  neigh_vb(NULL),

  dim_(0)

{
}


Element::Element(unsigned int dim, Mesh *mesh_in, Region reg)
{
    init(dim, mesh_in, reg);
}



void Element::init(unsigned int dim, Mesh *mesh_in, Region reg) {
    pid=0;
    material=NULL;
    n_neighs_vb=0;
    neigh_vb=NULL;
    dim_=dim;
    mesh_=mesh_in;
    region_=reg;

    node = new Node * [ n_nodes()];
    edge_idx_ = new unsigned int [ n_sides()];
    boundary_idx_ = NULL;

    FOR_ELEMENT_SIDES(this, si) {
        edge_idx_[ si ]=Mesh::undef_idx;
    }
}


/**
 * SET THE "VOLUME" FIELD IN STRUCT ELEMENT
 */
double Element::volume() {
    double volume = measure() * material->size;
    //INPUT_CHECK(!(volume < NUM_ZERO),
    //        "Volume of the element is nearly zero (volume= %g)\n", volume);
    return volume;
}

/**
 * SET THE "METRICS" FIELD IN STRUCT ELEMENT
 */
double Element::measure() {
    switch (dim()) {
        case 0:
            return 1.0;
            break;
        case 1:
            return arma::norm(*(node[ 1 ]) - *(node[ 0 ]) , 2);
            break;
        case 2:
            return
                arma::norm(
                    arma::cross(*(node[1]) - *(node[0]), *(node[2]) - *(node[0])),
                    2
                ) / 2.0 ;
            break;
        case 3:
            return fabs(
                arma::dot(
                    arma::cross(*node[1] - *node[0], *node[2] - *node[0]),
                    *node[3] - *node[0] )
                ) / 6.0;
            break;
    }
}


/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */

arma::vec3 Element::centre() {
    int li;

    arma::vec3 centre;
    centre.zeros();

    FOR_ELEMENT_NODES(this, li) {
        centre += node[ li ]->point();
    }
    centre /= (double) n_nodes();
    //DBGMSG("%d: %f %f %f\n",ele.id(),ele->centre[0],ele->centre[1],ele->centre[2]);
    return centre;
}

/**
 * Count element sides of the space dimension @p side_dim.
 */
unsigned int Element::n_sides_by_dim(int side_dim)
{
    if (side_dim == dim()) return 1;

    unsigned int n = 0;
    for (unsigned int i=0; i<n_sides(); i++)
        if (side(i)->dim() == side_dim) n++;
    return n;
}

/**
 * Return pointer to @p nth side/node/element (depending on the dimension @p side_dim).
 */
/*
void *Element::side_by_dim(int side_dim, unsigned int n)
{
    if (side_dim == 0)
    {
        // TODO: Maybe here we should also return a side?

        return node[n];
    }
    else if (side_dim == dim)
    {
        ASSERT(n==0, "Number of side is out of range.");
        return this;
    }
    else
    {
        unsigned int count = 0;
        for (unsigned int i=0; i<n_sides(); i++)
        {
            if (side(i)->dim() == side_dim)
            {
                if (count == n)
                {
                    return side(i);
                }
                else
                {
                    count++;
                }
            }
        }
        xprintf(Warn, "Side not found.");
    }
}
*/

/**
 * Return pointer to the @p node_id-th node of the side.
 */
const Node *Element::side_node(int side_dim, unsigned int side_id, unsigned node_id)
{
    if (side_dim == 0)
    {
        return node[side_id];
    }
    else if (side_dim == dim())
    {
        ASSERT(side_id==0, "Number of side is out of range.");
        return this->node[node_id];
    }
    else
    {
        unsigned int count = 0;
        for (unsigned int i=0; i<n_sides(); i++)
        {
            if (side(i)->dim() == side_dim)
            {
                if (count == side_id)
                {
                    return side(i)->node(node_id);
                }
                else
                {
                    count++;
                }
            }
        }
        xprintf(Warn, "Side not found.");
    }
}

ElementAccessor< 3 > Element::element_accessor()
{
  return mesh_->element_accessor( mesh_->element.index(this) );
}



#if 0

/**
 * make_block_d(ElementFullIter ele)
 */
void make_block_d(Mesh *mesh, ElementFullIter ele) {
    F_ENTRY;

    int ngi, iCol;
    struct Neighbour *ngh;
    ElementFullIter ele2 = ELEMENT_FULL_ITER_NULL(mesh);

    ele->d_row_count = 1 ;//+ ele->n_neighs_vv; // diagonal allways + noncompatible neighbours
    ele->d_col = (int*) xmalloc(ele->d_row_count * sizeof ( int));
    ele->d_el = (int*) xmalloc(ele->d_row_count * sizeof ( int));
    ele->d_val = (double*) xmalloc(ele->d_row_count * sizeof ( double));

    // set diagonal	on zero positon (D_DIAG == 0)
    ele->d_col[D_DIAG] = ele->b_row;
    ele->d_el[D_DIAG] = ele.index();
    ele->d_val[D_DIAG] = 0.0;

    // "Compatible" neighbours of higher dimensions

    FOR_ELM_NEIGHS_VB(ele, ngi) {
        ngh = ele->neigh_vb[ ngi ];
        ele->d_val[ D_DIAG ] -= ngh->sigma * ngh->side[1]->metric();
    }
    iCol = 1;

    // "Noncompatible" neighbours
/*
    FOR_ELM_NEIGHS_VV(ele, ngi) {
        ngh = ele->neigh_vv[ ngi ];
        // get neigbour element, and set appropriate column
        DBGMSG(" el1: %p el0: %p",ngh->element[1], ngh->element[0]);
        ele2 = ELEMENT_FULL_ITER(mesh,  (ngh->element[ 0 ] == ele) ? ngh->element[ 1 ] : ngh->element[ 0 ] );
        ele->d_el[ iCol ] = ele2.index();
        ele->d_col[ iCol ] = ele2->b_row;

        // add both sides of comunication
        double measure;
        if (ele->dim < ele2->dim) {
            measure = ele->measure;
        } else {
            measure = ele2->measure;
        }
        //DBGMSG("meas: %g\n",measure );
        ele->d_val[ D_DIAG ] -= ngh->sigma * ngh->geom_factor*measure;
        ele->d_val[ iCol ] += ngh->sigma * ngh->geom_factor*measure;

        iCol++;
    }*/
}

/**
 * make_block_e(ElementFullIter ele)
 */
/*
void make_block_e(ElementFullIter ele, Mesh *mesh) {
    int ngi, ci;
    struct Neighbour *ngh;

    ele->e_row_count = ele->n_neighs_vb;
    if (ele->e_row_count == 0) return;
    // alloc
    ele->e_col = (int*) xmalloc(ele->e_row_count * sizeof ( int));
    ele->e_edge_idx = (int*) xmalloc(ele->e_row_count * sizeof ( int));
    ele->e_val = (double*) xmalloc(ele->e_row_count * sizeof ( double));

    ci = 0;

    FOR_ELM_NEIGHS_VB(ele, ngi) {
        ngh = ele->neigh_vb[ ngi ];
        ele->e_col[ ci ] = ngh->edge->c_row;
        ele->e_val[ ci ] = ngh->sigma * ngh->side[1]->metric(); //DOPLNENO   * ngh->side[1]->metrics
        ele->e_edge_idx[ci] = mesh->edge.index(ngh->edge);
        ci++;
    }
}*/

/**
 * gets max,min, abs max, abs min of all local matrices
 * NEVER USED - may not work
 */
void block_A_stats(Mesh* mesh) {

    int i;
    double *loc;
    double a_min, a_max;
    double a_abs_min, a_abs_max;

    a_min = 1e32;
    a_max = -1e32;
    a_abs_min = 1e32;
    a_abs_max = -1e32;

    FOR_ELEMENTS(mesh, ele) {
        for (loc = ele->loc, i = 0; i < ele->n_sides * ele->n_sides; i++) {
            if (loc[i] < a_min) a_min = loc[i];
            if (loc[i] > a_max) a_max = loc[i];
            if (fabs(loc[i]) < a_abs_min) a_abs_min = fabs(loc[i]);
            if (fabs(loc[i]) > a_abs_max) a_abs_max = fabs(loc[i]);
            if (fabs(loc[i]) > 1e3)
                xprintf(Msg, "Big number: eid:%d area %g\n", ele.id(), ele->measure());
        }
    }

    xprintf(MsgVerb, "Statistics of the block A:\n")/*orig verb 6*/;
    xprintf(MsgVerb, "Minimal value: %g\tMaximal value: %g\n", a_min, a_max)/*orig verb 6*/;
    xprintf(MsgVerb, "Minimal absolute value: %g\tMaximal absolute value: %g\n", a_abs_min, a_abs_max)/*orig verb 6*/;
}

/**
 * gets max,min, abs max, abs min of all diagonals of local matrices
 * NEVER USED - may not work
 */
void diag_A_stats(Mesh* mesh) {

    int i;
    double *loc;
    double a_min, a_max;
    double a_abs_min, a_abs_max;

    a_min = 1e32;
    a_max = -1e32;
    a_abs_min = 1e32;
    a_abs_max = -1e32;

    FOR_ELEMENTS(mesh, ele) {
        for (loc = ele->loc, i = 0; i < ele->n_sides * ele->n_sides; i += ele->n_sides + 1) {
            // go through diagonal of ele->loc
            if (loc[i] < a_min) a_min = loc[i];
            if (loc[i] > a_max) a_max = loc[i];
            if (fabs(loc[i]) < a_abs_min) a_abs_min = fabs(loc[i]);
            if (fabs(loc[i]) > a_abs_max) a_abs_max = fabs(loc[i]);
        }
    }

    xprintf(MsgVerb, "Statistics of the diagonal of the block A:\n")/*orig verb 6*/;
    xprintf(MsgVerb, "Minimal value: %g\tMaximal value: %g\n", a_min, a_max)/*orig verb 6*/;
    xprintf(MsgVerb, "Minimal absolute value: %g\tMaximal absolute value: %g\n", a_abs_min, a_abs_max)/*orig verb 6*/;
}
#endif

//-----------------------------------------------------------------------------
// vim: set cindent:
