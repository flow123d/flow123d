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
 *
 * @file
 * @ingroup mesh
 * @brief Some side related functions - should be made strictly geometric.
 *
 */

#include "system/system.hh"
#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "sides.h"

// following deps. should be removed
#include "problem.h"
#include "mesh/boundaries.h"
#include "materials.hh"
#include "transport.h"
#include "constantdb.h"

static struct Side *new_side(void);
static void add_to_side_list(Mesh*, struct Side*);
static void init_side(struct Side*);
static int count_sides(Mesh*);
static void calc_side_c_row(struct Side*);
static void calc_side_c_col(Mesh*);
static void calc_side_c_val(struct Side*);
static void calc_side_rhs(struct Side*);
//static void calc_side_rhs_dens(struct Side*, struct Problem*, Mesh*);
static double side_length_line(struct Side*);
static double side_area_triangle(struct Side*);
static void calc_side_normal(struct Side*);
static void side_normal_point(struct Side*);
static void side_normal_line(struct Side*);
static void side_normal_triangle(struct Side*);
static void calc_side_centre(struct Side*);
static void side_centre_point(struct Side*);
static void side_centre_line(struct Side*);
static void side_centre_triangle(struct Side*);

//=============================================================================
// CREATE AND PREFILL LIST OF SIDES
//=============================================================================

void make_side_list(Mesh* mesh) {
    F_ENTRY;

    int si;
    struct Side *sde;

    ASSERT(!(mesh == NULL), "NULL as argument of function make_side_list()\n");
    xprintf(Msg, "Creating sides...")/*orig verb 2*/;
    mesh->n_sides = count_sides(mesh);
    for (si = 0; si < mesh->n_sides; si++) {
        sde = new_side();
        ASSERT(!(sde == NULL), "Cannot create side %d\n", si);
        add_to_side_list(mesh, sde);
        sde->id = si;
    }
    xprintf(MsgVerb, " O.K. %d sides created.", mesh->n_sides)/*orig verb 4*/;
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
//
//=============================================================================

int count_sides(Mesh* mesh) {
    F_ENTRY;

    int rc = 0;

    FOR_ELEMENTS(ele) {
        rc += ele->n_sides;
    }
    return rc;
}
//=============================================================================
// CREATE NEW SIDE
//=============================================================================

struct Side *new_side(void) {
    struct Side *sde;

    sde = (struct Side*) xmalloc(sizeof ( struct Side));
    init_side(sde);
    return sde;
}
//=============================================================================
// INIT DATA OF PARTICULAR SIDE
//=============================================================================

void init_side(struct Side *sde) {
    ASSERT(!(sde == NULL), "NULL as argument of function init_side()\n");
    sde->id = NDEF;
    sde->type = NDEF;
    sde->shape = NDEF;
    sde->dim = NDEF;
    sde->element = NULL;
    sde->lnum = NDEF;
    sde->n_nodes = NDEF;
    sde->node = NULL;
    sde->cond = NULL;
    sde->edge = NULL;
    sde->prev = NULL;
    sde->next = NULL;
    sde->neigh_bv = NULL;
    sde->metrics = 0.0;
    sde->normal[ 0 ] = 0.0;
    sde->normal[ 1 ] = 0.0;
    sde->normal[ 2 ] = 0.0;
    sde->centre[ 0 ] = 0.0;
    sde->centre[ 1 ] = 0.0;
    sde->centre[ 2 ] = 0.0;
    sde->c_row = NDEF;
    sde->c_col = NDEF;
    sde->c_val = 0.0;
    sde->flux = 0.0;
    sde->scalar = 0.0;
    sde->aux = NDEF;
    sde->faux = 0.0;
}
//=============================================================================
//
//=============================================================================

void add_to_side_list(Mesh* mesh, struct Side* sde) {
    F_ENTRY;

    ASSERT(!((mesh == NULL) || (sde == NULL)), "NULL as an argument of function add_to_side_list()\n");
    // First side in the list
    if (mesh->side == NULL && mesh->l_side == NULL) {
        mesh->side = sde;
        mesh->l_side = sde;
        sde->prev = NULL;
        sde->next = NULL;
        return;
    }
    // If something is wrong with the list
    ASSERT(!((mesh->side == NULL) || (mesh->l_side == NULL)), "Inconsistency in the side list\n");
    // Add after last side
    sde->next = NULL;
    sde->prev = mesh->l_side;
    mesh->l_side->next = sde;
    mesh->l_side = sde;
}
//=============================================================================
// CALCULATE PROPERTIES OF ALL SIDES OF THE MESH
//=============================================================================

void side_calculation_mh(Mesh* mesh) {
    struct Side *sde;

    xprintf(Msg, "Calculating properties of sides... ")/*orig verb 2*/;
    ASSERT(NONULL(mesh), "NULL as 'mesh' argument.\n");
    ASSERT(mesh->n_sides != NDEF && NONULL(mesh->side), "No side list.\n");
    calc_side_c_col(mesh);

    FOR_SIDES(sde) {
        calc_side_c_row(sde);
        calc_side_c_val(sde);
        calc_side_metrics(sde);
        calc_side_normal(sde);
        calc_side_centre(sde);
/*
        if (ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY)
            calc_side_rhs_dens(sde, problem, mesh);
        else */
            calc_side_rhs(sde);
    }
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// FILL THE "C_ROW" FIELD IN STRUCT SIDE FOR SIDE
//=============================================================================

void calc_side_c_row(struct Side *sde) {
    sde->c_row = sde->edge->c_row;
}
//=============================================================================
// FILL THE "C_COL" FIELD IN STRUCT SIDE FOR SIDE
//=============================================================================

void calc_side_c_col(Mesh* mesh) {
    //ElementIter ele;
    int li, i;

    i = 0;
    FOR_ELEMENTS(ele)
    for (li = 0; li < ele->n_sides; li++)
        ele->side[ li ]->c_col = i++;
}
//=============================================================================
// CALCULATE METRICS OF THE SIDE
//=============================================================================

void calc_side_metrics(struct Side *sde) {
    switch (sde->shape) {
        case S_POINT:
            sde->metrics = 1.0 * sde->element->material->size; //UPDATE
            break;
        case S_LINE:
            sde->metrics = side_length_line(sde) * sde->element->material->size; //UPDATE
            break;
        case S_TRIANGLE:
            sde->metrics = side_area_triangle(sde);
            break;
    }
}
//=============================================================================
//
//=============================================================================

double side_length_line(struct Side *sde) {
    double rc, u[ 3 ];

    u[ 0 ] = sde->node[ 1 ]->getX() - sde->node[ 0 ]->getX();
    u[ 1 ] = sde->node[ 1 ]->getY() - sde->node[ 0 ]->getY();
    u[ 2 ] = sde->node[ 1 ]->getZ() - sde->node[ 0 ]->getZ();
    rc = sqrt(u[ 0 ] * u[ 0 ] + u[ 1 ] * u[ 1 ] + u[ 2 ] * u[ 2 ]);
    return rc;
}
//=============================================================================
//
//=============================================================================

double side_area_triangle(struct Side *sde) {
    double u[ 3 ], v[ 3 ], n[ 3 ];
    double rc;

    u[ 0 ] = sde->node[ 1 ]->getX() - sde->node[ 0 ]->getX();
    u[ 1 ] = sde->node[ 1 ]->getY() - sde->node[ 0 ]->getY();
    u[ 2 ] = sde->node[ 1 ]->getZ() - sde->node[ 0 ]->getZ();
    v[ 0 ] = sde->node[ 2 ]->getX() - sde->node[ 0 ]->getX();
    v[ 1 ] = sde->node[ 2 ]->getY() - sde->node[ 0 ]->getY();
    v[ 2 ] = sde->node[ 2 ]->getZ() - sde->node[ 0 ]->getZ();
    vector_product(u, v, n);
    rc = fabs(0.5 * sqrt(n[ 0 ] * n[ 0 ] + n[ 1 ] * n[ 1 ] +
            n[ 2 ] * n[ 2 ]));
    return rc;
}
//=============================================================================
// CALCULATE NORMAL OF THE SIDE
//=============================================================================

void calc_side_normal(struct Side *sde) {
    switch (sde->shape) {
        case S_POINT:
            side_normal_point(sde);
            break;
        case S_LINE:
            side_normal_line(sde);
            break;
        case S_TRIANGLE:
            side_normal_triangle(sde);
            break;
    }
}
//=============================================================================
//
//=============================================================================

void side_normal_point(struct Side *sde) {
    ElementIter ele;

    ele = sde->element;
    sde->normal[ 0 ] = ele->node[ 1 ]->getX() - ele->node[ 0 ]->getX();
    sde->normal[ 1 ] = ele->node[ 1 ]->getY() - ele->node[ 0 ]->getY();
    sde->normal[ 2 ] = ele->node[ 1 ]->getZ() - ele->node[ 0 ]->getZ();
    normalize_vector(sde->normal);
    if (sde->node[ 0 ] == ele->node[ 0 ])
        scale_vector(sde->normal, -1);
}
//=============================================================================
//
//=============================================================================

void side_normal_line(struct Side *sde) {
    ElementIter ele;
    double s[ 3 ];
    double in[ 3 ];
    double en[ 3 ], u[ 3 ], v[ 3 ];

    // At first, we need vector of the normal of the element
    ele = sde->element;
    u[ 0 ] = ele->node[ 1 ]->getX() - ele->node[ 0 ]->getX();
    u[ 1 ] = ele->node[ 1 ]->getY() - ele->node[ 0 ]->getY();
    u[ 2 ] = ele->node[ 1 ]->getZ() - ele->node[ 0 ]->getZ();
    v[ 0 ] = ele->node[ 2 ]->getX() - ele->node[ 0 ]->getX();
    v[ 1 ] = ele->node[ 2 ]->getY() - ele->node[ 0 ]->getY();
    v[ 2 ] = ele->node[ 2 ]->getZ() - ele->node[ 0 ]->getZ();
    vector_product(u, v, en);
    normalize_vector(en);
    // Now we can calculate the "normal" of our side
    s[ 0 ] = sde->node[ 1 ]->getX() - sde->node[ 0 ]->getX();
    s[ 1 ] = sde->node[ 1 ]->getY() - sde->node[ 0 ]->getY();
    s[ 2 ] = sde->node[ 1 ]->getZ() - sde->node[ 0 ]->getZ();
    vector_product(s, en, sde->normal);
    normalize_vector(sde->normal);
    in[ 0 ] = ele->centre[ 0 ] - sde->node[ 0 ]->getX();
    in[ 1 ] = ele->centre[ 1 ] - sde->node[ 0 ]->getY();
    in[ 2 ] = ele->centre[ 2 ] - sde->node[ 0 ]->getZ();
    if (scalar_product(sde->normal, in) > 0.0)
        scale_vector(sde->normal, -1.0);
}
//=============================================================================
//
//=============================================================================

void side_normal_triangle(struct Side *sde) {
    ElementIter ele;
    double u[ 3 ], v[ 3 ], in[ 3 ];

    ele = sde->element;
    u[ 0 ] = sde->node[ 1 ]->getX() - sde->node[ 0 ]->getX();
    u[ 1 ] = sde->node[ 1 ]->getY() - sde->node[ 0 ]->getY();
    u[ 2 ] = sde->node[ 1 ]->getZ() - sde->node[ 0 ]->getZ();
    v[ 0 ] = sde->node[ 2 ]->getX() - sde->node[ 0 ]->getX();
    v[ 1 ] = sde->node[ 2 ]->getY() - sde->node[ 0 ]->getY();
    v[ 2 ] = sde->node[ 2 ]->getZ() - sde->node[ 0 ]->getZ();
    vector_product(u, v, sde->normal);
    normalize_vector(sde->normal);
    in[ 0 ] = ele->centre[ 0 ] - sde->node[ 0 ]->getX();
    in[ 1 ] = ele->centre[ 1 ] - sde->node[ 0 ]->getY();
    in[ 2 ] = ele->centre[ 2 ] - sde->node[ 0 ]->getZ();
    if (scalar_product(sde->normal, in) > 0.0)
        scale_vector(sde->normal, -1.0);
}
//=============================================================================
// CALCULATE VALUE IN THE BLOCK C
//=============================================================================

void calc_side_c_val(struct Side *sde) {
    sde->c_val = 1.0;
    if (sde->cond == NULL)
        return;
    if (sde->cond->type == DIRICHLET)
        sde->c_val = 0.0;
}
//=============================================================================


//=============================================================================
// CALCULATE VALUE ON THE RHS
//=============================================================================

void calc_side_rhs(struct Side *sde) {
    ElementIter ele;

    ASSERT(!(sde == NULL), "NULL argument\n", sde->id);
    ele = sde->element;
    ASSERT(!(ele == NULL), "Element of the side %d not defined\n", sde->id);
    ele->rhs[ sde->lnum ] += (ele->centre[ 2 ] - sde->centre[ 2 ]);
}
//=============================================================================

//======BP F.��r===============================================================
// CALCULATE VALUE ON THE RHS (Density)
//=============================================================================
void calc_side_rhs_dens(struct Side* sde, struct Problem* problem, Mesh* mesh) {
    /*

    ASSERT(!((sde == NULL) || (problem == NULL) || (mesh == NULL)), "NULL argument to calc_side_rhs_dens()");

    Transport *transport = problem->transport;
    int n_subst = mesh->n_substances;
    ElementFullIter ele = ELEMENT_FULL_ITER(sde->element);

    ASSERT(!(ele == NULL), "Element of the side %d not defined\n", sde->id);

    // compute total density of dissolved matter
    double sss = 0.0;
    for (int sbi = 0; sbi < n_subst; sbi++)
        sss += transport->substance_density_scale[sbi] * transport->out_conc[MOBILE][sbi][ele.index()];

    //xprintf( MsgVerb, " %f ", ele->start_conc->conc[0] );
    ele->rhs[sde->lnum] += (ele->centre[2] - sde->centre[2]) * (1 + sss
            / (ConstantDB::getInstance()->getDouble("Rho")));
    // * (1 + ele->start_conc->conc[0] / ConstantDB::getInstance()->getDouble("Rho"));
    //xprintf( MsgVerb, " %f ", (ele->centre[ 2 ] - sde->centre[ 2 ])
    /*  for( sbi = 0; sbi < n_subst; sbi++ ) {
     xprintf( MsgVerb, "%d %f %f \n", sbi, ele->conc [0] , ele->start_conc->conc[ sbi ] );
     }        */

}
//=============================================================================

//=============================================================================

void side_shape_specific(Mesh* mesh) {
    struct Side *sde;
    ElementIter ele;

    xprintf(MsgVerb, "   Filling shape specific data for sides... ")/*orig verb 5*/;
    ASSERT(NONULL(mesh), "NULL as mesh argument!\n");
    ASSERT((mesh->n_sides != NDEF) && NONULL(mesh->side), "No side list.\n");

    FOR_SIDES(sde) {
        ele = sde->element;
        ASSERT(!(ele == NULL), "Side %d has no reference to its element\n", sde->id);
        switch (ele->type) {
            case LINE:
                sde->shape = S_POINT;
                sde->dim = 0;
                sde->n_nodes = 1;
                break;
            case TRIANGLE:
                sde->shape = S_LINE;
                sde->dim = 1;
                sde->n_nodes = 2;
                break;
            case TETRAHEDRON:
                sde->shape = S_TRIANGLE;
                sde->dim = 2;
                sde->n_nodes = 3;
                break;
        }
        sde->node = (Node**) xmalloc(sde->n_nodes * sizeof (Node*));
    }
    xprintf(MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
// CALCULATE CENTRE OF THE SIDE
//=============================================================================

void calc_side_centre(struct Side *sde) {
    switch (sde->shape) {
        case S_POINT:
            side_centre_point(sde);
            break;
        case S_LINE:
            side_centre_line(sde);
            break;
        case S_TRIANGLE:
            side_centre_triangle(sde);
            break;
    }
}
//=============================================================================
//
//=============================================================================

void side_centre_point(struct Side *sde) {
    sde->centre[ 0 ] = sde->node[ 0 ]->getX();
    sde->centre[ 1 ] = sde->node[ 0 ]->getY();
    sde->centre[ 2 ] = sde->node[ 0 ]->getZ();
}
//=============================================================================
//
//=============================================================================

void side_centre_line(struct Side *sde) {
    sde->centre[ 0 ] = (sde->node[ 0 ]->getX() + sde->node[ 1 ]->getX()) / 2.0;
    sde->centre[ 1 ] = (sde->node[ 0 ]->getY() + sde->node[ 1 ]->getY()) / 2.0;
    sde->centre[ 2 ] = (sde->node[ 0 ]->getZ() + sde->node[ 1 ]->getZ()) / 2.0;
}
//=============================================================================
//
//=============================================================================

void side_centre_triangle(struct Side *sde) {
    sde->centre[ 0 ] = (sde->node[ 0 ]->getX() +
            sde->node[ 1 ]->getX() +
            sde->node[ 2 ]->getX()) / 3.0;
    sde->centre[ 1 ] = (sde->node[ 0 ]->getY() +
            sde->node[ 1 ]->getY() +
            sde->node[ 2 ]->getY()) / 3.0;
    sde->centre[ 2 ] = (sde->node[ 0 ]->getZ() +
            sde->node[ 1 ]->getZ() +
            sde->node[ 2 ]->getZ()) / 3.0;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
