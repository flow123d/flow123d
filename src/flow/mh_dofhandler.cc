/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    mh_dofhandler.cc
 * @brief   
 */

#include "flow/mh_dofhandler.hh"
#include "mesh/mesh.h"
#include "mesh/side_impl.hh"

void MH_DofHandler::reinit(Mesh *mesh) {
    elem_side_to_global.resize(mesh->n_elements() );
    FOR_ELEMENTS(mesh, ele) elem_side_to_global[ele.index()].resize(ele->n_sides());

    unsigned int i_side_global=0;
    FOR_ELEMENTS(mesh, ele) {
        for(unsigned int i_lside=0; i_lside < ele->n_sides(); i_lside++)
            elem_side_to_global[ele.index()][i_lside] = i_side_global++;
    }
}


unsigned int MH_DofHandler::side_dof(const SideIter side) const {
    return elem_side_to_global[ side->element().index() ][ side->el_idx() ];
}


void MH_DofHandler::set_solution( double time, double * solution, double precision) {
    ASSERT( solution != NULL, "Empty solution.\n");
    mh_solution = solution;
    solution_precision = precision;
    time_ = time;
}

/// temporary replacement for DofHandler accessor, flux through given side
double MH_DofHandler::side_flux(const Side &side) const {
    return mh_solution[ elem_side_to_global[ side.element().index() ][ side.el_idx() ] ];
}

/// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
double MH_DofHandler::side_scalar(const Side &side) const {
    unsigned int i_edg = side.edge_idx();
    return mh_solution[ side.mesh()->n_sides() + side.mesh()->n_elements() + i_edg ];
}


double MH_DofHandler::element_scalar( ElementFullIter &ele ) const {
    return mh_solution[ ele->mesh_->n_sides() + ele.index() ];
}
