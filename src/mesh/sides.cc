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




//static void calc_side_rhs_dens(struct Side*, struct Problem*, Mesh*);







//=============================================================================
// CALCULATE PROPERTIES OF ALL SIDES OF THE MESH
//=============================================================================


//=============================================================================
// CALCULATE METRICS OF THE SIDE
//=============================================================================

double Side::measure() const {
    switch ( dim() ) {
        case 0:
            return 1.0;
        case 1: {
            arma::vec3 diff = node(1)->point();
            diff -= node(0)->point();
            return arma::norm( diff , 2 );
        }
        case 2: {
            arma::vec3 diff0 = node(1)->point() - node(0)->point();
            arma::vec3 diff1 = node(2)->point() - node(0)->point();
            return 0.5*arma::norm( arma::cross(diff0, diff1), 2);
        }
    }
}

//=============================================================================
// CALCULATE NORMAL OF THE SIDE
//=============================================================================

arma::vec3 Side::normal() const {
    switch ( dim() ) {
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

arma::vec3 Side::normal_point() const {
    ElementIter ele = element_;

    arma::vec3 normal(ele->node[1]->point());
    normal -= ele->node[0] ->point();

    normal /=arma::norm(normal,2);
    if ( node( 0 ) == ele->node[ 0 ] )
        return -normal;
    else
        return normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_line() const {
    ElementIter ele=element_;

    // At first, we need vector of the normal of the element
    arma::vec3 elem_normal=arma::cross( ele->node[1]->point() - ele->node[0]->point(),
                                        ele->node[2]->point() - ele->node[0]->point() );
    elem_normal /= norm( elem_normal, 2);

    // Now we can calculate the "normal" of our side
    arma::vec3 side_normal = arma::cross( node(1)->point() - node(0)->point() , elem_normal );
    side_normal /= norm( side_normal, 2);

    if ( dot( side_normal, ele->centre() - node(0)->point() ) > 0.0)
        return -side_normal;
    else
        return side_normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_triangle() const {
    ElementIter ele=element_;
    double u[ 3 ], v[ 3 ], in[ 3 ], normal[3];

    arma::vec3 side_normal=arma::cross( node(1)->point() - node(0)->point(),
                                        node(2)->point() - node(0)->point() );
    side_normal /= norm( side_normal, 2);

    in[ 0 ] = ele->centre()[ 0 ] - node( 0 )->getX();
    in[ 1 ] = ele->centre()[ 1 ] - node( 0 )->getY();
    in[ 2 ] = ele->centre()[ 2 ] - node( 0 )->getZ();
    if ( dot(side_normal, ele->centre() - node(0)->point() ) > 0.0)
        return -side_normal;
    else
        return side_normal;
}

//=============================================================================
// CALCULATE CENTRE OF THE SIDE
//=============================================================================

arma::vec3 Side::centre() const {
    arma::vec3 barycenter;
    barycenter.zeros();

    for(unsigned int i=0; i < n_nodes() ; i++)
        barycenter += node( i )->point();

    barycenter /= (double) n_nodes();
    return barycenter;
}


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

//-----------------------------------------------------------------------------
// vim: set cindent:
