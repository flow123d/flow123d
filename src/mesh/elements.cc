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

#include <strings.h>

#include "system/system.hh"
#include "xio.h"
#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "elements.h"

// following deps. should be removed
#include "problem.h"
#include "mesh/boundaries.h"
#include "materials.hh"

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"


static void calc_a_row(Mesh*);
static void calc_b_row(Mesh*);
//static ElementIter new_element(void);
//static void add_to_element_list(Mesh*, ElementIter);
static void make_block_e(ElementFullIter, Mesh *mesh );
//static void alloc_and_init_block_e(ElementIter );
static char supported_element_type(int);
static void element_type_specific(ElementFullIter );
static void element_allocation_independent(ElementFullIter );
static void make_block_d(Mesh *mesh, ElementFullIter );
static void calc_metrics(ElementFullIter );
static void calc_volume(ElementFullIter );
static double element_length_line(ElementFullIter );
static double element_area_triangle(ElementFullIter );
static double element_volume_tetrahedron(ElementFullIter );
static void calc_centre(ElementFullIter );
static void calc_rhs(ElementFullIter );
static void calc_rhs_b(ElementFullIter );
static void dirichlet_elm(ElementFullIter );

static void parse_element_properties_line(char*);
static void block_A_stats(Mesh*);
static void diag_A_stats(Mesh*);

Element::Element()
: type(0),
  mid(0),
  rid(0),
  pid(0),

  dim(0),
  n_sides(0),

  n_nodes(0),
  node(NULL),

  material(NULL),
  side(NULL),
  n_neighs_vv(0),
  neigh_vv(NULL),
  n_neighs_vb(0),
  neigh_vb(NULL),
 //*start_conc,
 //n_subst,
            // Material properties
  k(NULL),
  a(NULL),
  //stor(0),
            // Geometrical properties
  measure(0),
  volume(0),
            // Parameters of the basis functions
 bas_alfa(NULL),
 bas_beta(NULL),
 bas_gama(NULL),
 bas_delta(NULL),
            // Matrix
 loc(NULL),
 loc_inv(NULL),
 rhs(NULL),

 a_row(0),
 b_row(0),
 d_row_count(0),
 d_col(NULL),
 d_val(NULL),
 d_el(NULL),

 e_row_count(0),
 e_col(NULL),
 e_edge_idx(NULL),
 e_val(NULL),

 v_length(0),
 scalar(0),
 pscalar(0),
 balance(0),

 scalar_it(0),
// *conc_prev,
// *conc_prev_immobile,
// *conc_prev_sorb,
// *conc_prev_immobile_sorb,

 aux(0),
 faux(0)

{
 centre.zeros();
}

/**
 * CALCULATE PROPERTIES OF ALL ELEMENTS OF THE MESH
 */
void element_calculation_mh(Mesh* mesh) {

    F_ENTRY;

    ASSERT(NONULL(mesh), "No mesh for problem\n");
    ASSERT(mesh->element.size() > 0, "Empty mesh.\n");

    xprintf(Msg, "Calculating properties of elements... ")/*orig verb 2*/;

    calc_a_row(mesh);
    calc_b_row(mesh);

    FOR_ELEMENTS(ele) {
        if (ele->material->dimension != ele->dim) {
            xprintf(Warn, "Dimension %d of material doesn't match dimension %d of element %d.\n",
                    ele->material->dimension, ele->dim, ele.id());
        }
        ele->a = ele->material->hydrodynamic_resistence;

        calc_rhs(ele);
        dirichlet_elm(ele);
        make_block_d(mesh, ele);
        make_block_e(ele, mesh);
    }
    //block_A_stats( mesh );
    //diag_A_stats( mesh );
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}

/**
 * CALCULATE PROPERTIES OF ALL ELEMENTS OF THE MESH
 */
void make_element_geometry() {
    xprintf(Msg, "Calculating properties of elements... ")/*orig verb 2*/;

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);


    ASSERT(NONULL(mesh), "No mesh for problem\n");
    ASSERT(mesh->element.size() > 0, "Empty mesh.\n");

    FOR_ELEMENTS(ele) {
        //DBGMSG("\n ele: %d \n",ele.id());
        //FOR_ELEMENTS(ele1) {
        //    printf("%d(%d) ",ele1.id(),ele1->type);
        //    ele1->bas_alfa[0]=1.0;
       // }
        calc_metrics(ele);
        calc_volume(ele);
        calc_centre(ele);
    }

    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}

/**
 * CALCULATE THE "A_ROW" FIELD IN STRUCT ELEMENT
 */
void calc_a_row(Mesh* mesh) {
    int last = 0;

    FOR_ELEMENTS(ele) {
        ele->a_row = last;
        last += ele->n_sides;
    }
}

/**
 * CALCULATE THE "B_ROW" FIELD IN STRUCT ELEMENT
 */
void calc_b_row(Mesh* mesh) {
    int last;

    last = mesh->n_sides;

    FOR_ELEMENTS(ele) {
        ele->b_row = last++;
    }
}

/**
 * SET THE "VOLUME" FIELD IN STRUCT ELEMENT
 */
void calc_volume(ElementFullIter ele) {
    ele->volume = ele->measure * ele->material->size; //UPDATE
    //        ele->volume = ele->metrics * ele->size;	   // JB version
    INPUT_CHECK(!(ele->volume < NUM_ZERO),
            "Volume of the element %d is nearly zero (volume= %g)\n", ele.id(), ele->volume);
}

/**
 * SET THE "METRICS" FIELD IN STRUCT ELEMENT
 */
void calc_metrics(ElementFullIter ele) {
    switch (ele->type) {
        case LINE:
            ele->measure = element_length_line(ele);
            break;
        case TRIANGLE:
            ele->measure = element_area_triangle(ele);
            break;
        case TETRAHEDRON:
            ele->measure = element_volume_tetrahedron(ele);
            break;
    }
}

/**
 * CALCULATE LENGTH OF LINEAR ELEMENT
 */

double element_length_line(ElementFullIter ele) {

    return arma::norm(*ele->node[ 1 ] - *ele->node[ 0 ] , 2);
}

/**
 * CALCULATE AREA OF TRIANGULAR ELEMENT
 */
double element_area_triangle(ElementFullIter ele) {

    return
        arma::norm(
            arma::cross(*ele->node[1] - *ele->node[0], *ele->node[2] - *ele->node[0]),
            2
        ) / 2.0 ;
}

/**
 * CALCULATE VOLUME OF TETRAHEDRA ELEMENT
 */
double element_volume_tetrahedron(ElementFullIter ele) {
    return fabs(
        arma::dot(
            arma::cross(*ele->node[1] - *ele->node[0], *ele->node[2] - *ele->node[0]),
            *ele->node[3] - *ele->node[0] )
        ) / 6.0;
}

/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */

void calc_centre(ElementFullIter ele) {
    int li;

    ele->centre.zeros();

    FOR_ELEMENT_NODES(ele, li) {
        ele->centre += ele->node[ li ]->point();
    }
    ele->centre /= (double) ele->n_nodes;
    //DBGMSG("%d: %f %f %f\n",ele.id(),ele->centre[0],ele->centre[1],ele->centre[2]);

/*
    ele->centre[ 0 ] = 0.0;
    ele->centre[ 1 ] = 0.0;
    ele->centre[ 2 ] = 0.0;

    FOR_ELEMENT_NODES(ele, li) {
        ele->centre[ 0 ] += ele->node[ li ]->getX();
        ele->centre[ 1 ] += ele->node[ li ]->getY();
        ele->centre[ 2 ] += ele->node[ li ]->getZ();
    }
    ele->centre[ 0 ] /= (double) ele->n_nodes;
    ele->centre[ 1 ] /= (double) ele->n_nodes;
    ele->centre[ 2 ] /= (double) ele->n_nodes;
*/
}


/**
 * SET THE "RHS[]" FIELD IN STRUCT ELEMENT
 */
void calc_rhs(ElementFullIter ele) {
    int li;

    FOR_ELEMENT_SIDES(ele, li) {
        ele->rhs[ li ] = 0.0;
    }
}

/**
 * CORRECT RHS IN CASE, WHEN DIRICHLET'S CONDITION IS GIVEN
 */
void dirichlet_elm(ElementFullIter ele) {
    int li;
    struct Boundary *bcd;

    FOR_ELEMENT_SIDES(ele, li) {
        bcd = ele->side[ li ]->cond;
        if (bcd == NULL) continue;
        if (bcd->type == DIRICHLET) ele->rhs[ li ] -= bcd->scalar;
    }
}


/**
 * make_block_d(ElementFullIter ele)
 */
void make_block_d(Mesh *mesh, ElementFullIter ele) {
    F_ENTRY;

    int ngi, iCol;
    struct Neighbour *ngh;
    ElementFullIter ele2 = ELEMENT_FULL_ITER_NULL;

    ele->d_row_count = 1 + ele->n_neighs_vv; // diagonal allways + noncompatible neighbours
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
        ele->d_val[ D_DIAG ] -= ngh->sigma * ngh->side[1]->metrics;
    }
    iCol = 1;

    // "Noncompatible" neighbours

    FOR_ELM_NEIGHS_VV(ele, ngi) {
        ngh = ele->neigh_vv[ ngi ];
        // get neigbour element, and set appropriate column
        DBGMSG(" el1: %p el0: %p",ngh->element[1], ngh->element[0]);
        ele2 = ELEMENT_FULL_ITER( (ngh->element[ 0 ] == ele) ? ngh->element[ 1 ] : ngh->element[ 0 ] );
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
    }
}

/**
 * make_block_e(ElementFullIter ele)
 */
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
        ele->e_val[ ci ] = ngh->sigma * ngh->side[1]->metrics; //DOPLNENO   * ngh->side[1]->metrics
        ele->e_edge_idx[ci] = mesh->edge.index(ngh->edge);
        ci++;
    }
}

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

    FOR_ELEMENTS(ele) {
        for (loc = ele->loc, i = 0; i < ele->n_sides * ele->n_sides; i++) {
            if (loc[i] < a_min) a_min = loc[i];
            if (loc[i] > a_max) a_max = loc[i];
            if (fabs(loc[i]) < a_abs_min) a_abs_min = fabs(loc[i]);
            if (fabs(loc[i]) > a_abs_max) a_abs_max = fabs(loc[i]);
            if (fabs(loc[i]) > 1e3)
                xprintf(Msg, "Big number: eid:%d area %g\n", ele.id(), ele->measure);
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

    FOR_ELEMENTS(ele) {
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

//-----------------------------------------------------------------------------
// vim: set cindent:
