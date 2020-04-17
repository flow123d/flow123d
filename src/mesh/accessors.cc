/*!
 *
 * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    accessors.cc
 * @brief   Implementation of other functions of the mesh accessors.
 */

#include "mesh/accessors.hh"
#include "system/system.hh"
#include "mesh/mesh.h"


//=============================================================================
// CALCULATE METRICS OF THE SIDE
//=============================================================================

double Side::measure() const {
    switch ( dim() ) {
        case 0:
            return 1.0;
        case 1: {
            return arma::norm( *node(1) - *node(0) , 2 );
        }
        case 2: {
            arma::vec3 diff0 = *node(1) - *node(0);
            arma::vec3 diff1 = *node(2) - *node(0);
            return 0.5*arma::norm( arma::cross(diff0, diff1), 2);
        }
    }

    return 0.0;
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

    return arma::vec3("0 0 0");
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_point() const {
    ElementAccessor<3> ele = element();

    arma::vec3 normal(*ele.node(1));
    normal -= *ele.node(0);

    normal /=arma::norm(normal, 2);
    if ( node(0) == ele.node(0) )
        return -normal;
    else
        return normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_line() const {
	ElementAccessor<3> ele = element();

    // At first, we need vector of the normal of the element
    arma::vec3 elem_normal=arma::cross( *ele.node(1) - *ele.node(0),
                                        *ele.node(2) - *ele.node(0) );
    elem_normal /= norm( elem_normal, 2);

    // Now we can calculate the "normal" of our side
    arma::vec3 side_normal = arma::cross( *node(1) - *node(0) , elem_normal );
    side_normal /= norm( side_normal, 2);

    if ( arma::dot( side_normal, ele.centre() - *node(0) ) > 0.0)
        return -side_normal;
    else
        return side_normal;
}
//=============================================================================
//
//=============================================================================

arma::vec3 Side::normal_triangle() const {
    arma::vec3 side_normal=arma::cross( *node(1) - *node(0),
                                        *node(2) - *node(0) );
    side_normal /= arma::norm( side_normal, 2);

    if ( arma::dot(side_normal, element().centre() - *node(0) ) > 0.0)
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
        barycenter += *node( i );

    barycenter /= (double) n_nodes();
    return barycenter;
}


//=============================================================================
// CALCULATE THE SIDE DIAMETER
//=============================================================================

double Side::diameter() const {
	if (dim() == 0)
    {
        return 1.0;
    }
    else
    {
    	double h = 0.0;
        for (unsigned int i=0; i<n_nodes(); i++)
            for (unsigned int j=i+1; j<n_nodes(); j++)
                h = max(h, arma::norm(*node(i) - *node(j), 2) );
        return h;
    }
}