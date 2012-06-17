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
#include "mesh/mesh_types.hh"

// following deps. should be removed
#include "mesh/boundaries.h"
#include "materials.hh"
//#include "transport.h"


static void calc_side_c_row(Side&);
static void calc_side_c_col(Mesh*);
static void calc_side_c_val(Side&);
static void calc_side_rhs(Side&);
//static void calc_side_rhs_dens(struct Side*, struct Problem*, Mesh*);



Side::Side() {
    element = NULL;
    node = NULL;
    cond = NULL;
    edge = NULL;
    neigh_bv = NULL;
}

void Side::reinit(ElementIter ele, unsigned int set_dim, int set_id, int set_lnum) {
    element = ele;
    dim= set_dim;
    id = set_id;
    lnum = set_lnum;

    if (node != NULL) delete [] node;
    node = new Node*[n_nodes()];
}




//=============================================================================
// CALCULATE PROPERTIES OF ALL SIDES OF THE MESH
//=============================================================================

void side_calculation_mh(Mesh* mesh) {
    struct Side *sde;

    xprintf(Msg, "Calculating properties of sides... ")/*orig verb 2*/;
    ASSERT(NONULL(mesh), "NULL as 'mesh' argument.\n");
    ASSERT(mesh->n_sides(), "No side list.\n");
    calc_side_c_col(mesh);

    FOR_SIDES(mesh, sde) {
        calc_side_c_row(*sde);
        calc_side_c_val(*sde);

/*
        if (ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY)
            calc_side_rhs_dens(sde, problem, mesh);
        else */
            calc_side_rhs(*sde);
    }
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// FILL THE "C_ROW" FIELD IN STRUCT SIDE FOR SIDE
//=============================================================================

void calc_side_c_row(Side &sde) {
    sde.c_row = sde.edge->c_row;
}
//=============================================================================
// FILL THE "C_COL" FIELD IN STRUCT SIDE FOR SIDE
//=============================================================================

void calc_side_c_col(Mesh* mesh) {
    //ElementIter ele;
    int li, i;

    i = 0;
    FOR_ELEMENTS(mesh, ele)
    for (li = 0; li < ele->n_sides; li++)
        ele->side[ li ]->c_col = i++;
}
//=============================================================================
// CALCULATE METRICS OF THE SIDE
//=============================================================================

double Side::metric() {
    switch (dim) {
        case 0:
            return 1.0 * element->material->size; //UPDATE
            break;
        case 1:
            return length_line() * element->material->size; //UPDATE
            break;
        case 2:
            return  area_triangle();
            break;
    }
}
//=============================================================================
//
//=============================================================================

double Side::length_line() {
    double rc, u[ 3 ];

    u[ 0 ] = node[ 1 ]->getX() - node[ 0 ]->getX();
    u[ 1 ] = node[ 1 ]->getY() - node[ 0 ]->getY();
    u[ 2 ] = node[ 1 ]->getZ() - node[ 0 ]->getZ();
    rc = sqrt(u[ 0 ] * u[ 0 ] + u[ 1 ] * u[ 1 ] + u[ 2 ] * u[ 2 ]);
    return rc;
}
//=============================================================================
//
//=============================================================================

double Side::area_triangle() {
    double u[ 3 ], v[ 3 ], n[ 3 ];
    double rc;

    u[ 0 ] = node[ 1 ]->getX() - node[ 0 ]->getX();
    u[ 1 ] = node[ 1 ]->getY() - node[ 0 ]->getY();
    u[ 2 ] = node[ 1 ]->getZ() - node[ 0 ]->getZ();
    v[ 0 ] = node[ 2 ]->getX() - node[ 0 ]->getX();
    v[ 1 ] = node[ 2 ]->getY() - node[ 0 ]->getY();
    v[ 2 ] = node[ 2 ]->getZ() - node[ 0 ]->getZ();
    vector_product(u, v, n);
    rc = fabs(0.5 * sqrt(n[ 0 ] * n[ 0 ] + n[ 1 ] * n[ 1 ] +
            n[ 2 ] * n[ 2 ]));
    return rc;
}
//=============================================================================
// CALCULATE NORMAL OF THE SIDE
//=============================================================================

arma::vec3 Side::normal() {
    switch (dim) {
        case 0:
            return normal_point();
        case 1:
            return normal_line();
        case 2:
            return normal_triangle();
    }
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_point() {
    ElementIter ele = element;

    arma::vec3 normal(ele->node[1]->point());
    normal -= ele->node[0] ->point();

    normal /=arma::norm(normal,2);
    if ( node[ 0 ] == ele->node[ 0 ] )
        return -normal;
    else
        return normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_line() {
    ElementIter ele=element;

    // At first, we need vector of the normal of the element
    arma::vec3 elem_normal=arma::cross( ele->node[1]->point() - ele->node[0]->point(),
                                        ele->node[2]->point() - ele->node[0]->point() );
    elem_normal /= norm( elem_normal, 2);

    // Now we can calculate the "normal" of our side
    arma::vec3 side_normal = arma::cross( node[1]->point() - node[0]->point() , elem_normal );
    side_normal /= norm( side_normal, 2);

    if ( dot( side_normal, ele->centre() - node[0]->point() ) > 0.0)
        return -side_normal;
    else
        return side_normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_triangle() {
    ElementIter ele=element;
    double u[ 3 ], v[ 3 ], in[ 3 ], normal[3];

    arma::vec3 side_normal=arma::cross( node[1]->point() - node[0]->point(),
                                        node[2]->point() - node[0]->point() );
    side_normal /= norm( side_normal, 2);

    in[ 0 ] = ele->centre()[ 0 ] - node[ 0 ]->getX();
    in[ 1 ] = ele->centre()[ 1 ] - node[ 0 ]->getY();
    in[ 2 ] = ele->centre()[ 2 ] - node[ 0 ]->getZ();
    if ( dot(side_normal, ele->centre() - node[0]->point() ) > 0.0)
        return -side_normal;
    else
        return side_normal;
}
//=============================================================================
// CALCULATE VALUE IN THE BLOCK C
//=============================================================================

void calc_side_c_val(Side &sde) {
    sde.c_val = 1.0;
    if (sde.cond == NULL)
        return;
    if (sde.cond->type == DIRICHLET)
        sde.c_val = 0.0;
}
//=============================================================================


//=============================================================================
// CALCULATE VALUE ON THE RHS -
//=============================================================================

void calc_side_rhs(Side &sde) {
    ElementIter ele;

    ele = sde.element;
    ASSERT(!(ele == NULL), "Element of the side %d not defined\n", sde.id);
    ele->rhs[ sde.lnum ] += (ele->centre()[ 2 ] - sde.centre()[ 2 ]);
    /*
     * prbably zero order approximation of :
     * Int_{El} z div Psi - Int_{Side} z Psi .dot. normal
     */

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
/*
void side_shape_specific(Mesh* mesh) {
    struct Side *sde;
    ElementIter ele;

    xprintf(MsgVerb, "   Filling shape specific data for sides... ");
    ASSERT(NONULL(mesh), "NULL as mesh argument!\n");
    ASSERT((mesh->n_sides != NDEF) && NONULL(mesh->side), "No side list.\n");

    FOR_SIDES(mesh, sde) {
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
    }
    xprintf(MsgVerb, "O.K.\n");
}*/
//=============================================================================
// CALCULATE CENTRE OF THE SIDE
//=============================================================================

arma::vec3 Side::centre() {
    arma::vec3 barycenter;
    barycenter.zeros();

    for(unsigned int i=0; i < n_nodes() ; i++)
        barycenter += node[ i ]->point();

    barycenter /= (double) n_nodes();
    return barycenter;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
