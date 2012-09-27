/*
 * mh_dofhandler.cc
 *
 *  Created on: Jun 18, 2012
 *      Author: jb
 */

#include "flow/mh_dofhandler.hh"
#include "mesh/mesh.h"
#include "mesh/side_impl.hh"

void MH_DofHandler::reinit(Mesh *mesh) {
    elem_side_to_global.resize(mesh->n_elements() );
    FOR_ELEMENTS(mesh, ele) elem_side_to_global[ele.index()].resize(ele->n_sides());

    unsigned int i_side_global=0;
    FOR_ELEMENTS(mesh, ele) {
        for(int i_lside=0; i_lside < ele->n_sides(); i_lside++)
            elem_side_to_global[ele.index()][i_lside] = i_side_global++;
    }
}


unsigned int MH_DofHandler::side_dof(const SideIter side) const {
    return elem_side_to_global[ side->element().index() ][ side->el_idx() ];
}


void MH_DofHandler::set_solution( double * solution) {
    ASSERT( solution != NULL, "Empty solution.\n");
    mh_solution = solution;
}

/// temporary replacement for DofHandler accessor, flux through given side
double MH_DofHandler::side_flux(const Side &side) const {
    return mh_solution[ elem_side_to_global[ side.element().index() ][ side.el_idx() ] ];
}

/// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
double MH_DofHandler::side_scalar(const Side &side) const {
    Edge * edg = side.edge();
    return mh_solution[ side.mesh()->n_sides() + side.mesh()->n_elements() + side.mesh()->edge.index( edg )];
}


double MH_DofHandler::element_scalar( ElementFullIter &ele ) const {
    return mh_solution[ ele->mesh_->n_sides() + ele.index() ];
}
