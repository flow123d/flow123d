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
#include "la/local_to_global_map.hh"
#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/neighbours.h"
#include "system/index_types.hh"
#include "system/sys_profiler.hh"

MH_DofHandler::MH_DofHandler()
:
el_4_loc(nullptr),
row_4_el(nullptr),
side_id_4_loc(nullptr),
side_row_4_id(nullptr),
edge_4_loc(nullptr),
row_4_edge(nullptr),
edge_ds(nullptr),
el_ds(nullptr),
side_ds(nullptr)
{}

MH_DofHandler::~MH_DofHandler()
{
    delete edge_ds;
    //  delete el_ds;
    delete side_ds;

    //  xfree(el_4_loc);
    delete [] row_4_el;
    delete [] side_id_4_loc;
    delete [] side_row_4_id;
    delete [] edge_4_loc;
    delete []  row_4_edge;
}


void MH_DofHandler::reinit(Mesh *mesh) {
    mesh_ = mesh;
    elem_side_to_global.resize(mesh->n_elements() );
    for (auto ele : mesh->elements_range()) elem_side_to_global[ ele.idx() ].resize(ele->n_sides());

    unsigned int i_side_global=0;
    for (auto ele : mesh->elements_range()) {
        for(unsigned int i_lside=0; i_lside < ele->n_sides(); i_lside++)
            elem_side_to_global[ ele.idx() ][i_lside] = i_side_global++;
    }

    prepare_parallel();

}




// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void MH_DofHandler::prepare_parallel() {

    START_TIMER("prepare parallel");

    LongIdx *loc_part; // optimal (edge,el) partitioning (local chunk)
    LongIdx *id_4_old; // map from old idx to ids (edge,el)
    int loc_i;

    int e_idx;


    //ierr = MPI_Barrier(PETSC_COMM_WORLD);
    //OLD_ASSERT(ierr == 0, "Error in MPI_Barrier.");

    // row_4_el will be modified so we make a copy of the array from mesh
    row_4_el = new LongIdx[mesh_->n_elements()];
    std::copy(mesh_->get_row_4_el(), mesh_->get_row_4_el()+mesh_->n_elements(), row_4_el);
    el_4_loc = mesh_->get_el_4_loc();
    el_ds = mesh_->get_el_ds();

    //optimal element part; loc. els. id-> new el. numbering
    Distribution init_edge_ds(DistributionLocalized(), mesh_->n_edges(), PETSC_COMM_WORLD);
    // partitioning of edges, edge belongs to the proc of his first element
    // this is not optimal but simple
    loc_part = new LongIdx[init_edge_ds.lsize()];
    id_4_old = new LongIdx[mesh_->n_edges()];
    {
        loc_i = 0;
        for(auto edg: mesh_->edge_range())
        {
            unsigned int i_edg = edg.idx();
            // partition
            e_idx = edg.side(0)->element().idx();
            if (init_edge_ds.is_local(i_edg)) {
                // find (new) proc of the first element of the edge
                loc_part[loc_i++] = el_ds->get_proc(row_4_el[e_idx]);
            }
            // id array
            id_4_old[i_edg] = i_edg;
        }
    }
    Partitioning::id_maps(mesh_->n_edges(), id_4_old, init_edge_ds, loc_part, edge_ds, edge_4_loc, row_4_edge);
    delete[] loc_part;
    delete[] id_4_old;

    // create map from mesh global edge id to new local edge id
    unsigned int loc_edge_idx=0;
    for (unsigned int i_el_loc = 0; i_el_loc < el_ds->lsize(); i_el_loc++) {
        auto ele = mesh_->element_accessor( el_4_loc[i_el_loc] );
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            unsigned int mesh_edge_idx= ele.side(i)->edge_idx();
            if ( edge_new_local_4_mesh_idx_.count(mesh_edge_idx) == 0 )
                // new local edge
                edge_new_local_4_mesh_idx_[mesh_edge_idx] = loc_edge_idx++;
        }
    }

    //optimal side part; loc. sides; id-> new side numbering
    Distribution init_side_ds(DistributionBlock(), mesh_->n_sides(), PETSC_COMM_WORLD);
    // partitioning of sides follows elements
    loc_part = new LongIdx[init_side_ds.lsize()];
    id_4_old = new LongIdx[mesh_->n_sides()];
    {
        int is = 0;
        loc_i = 0;
    	for (auto ele : mesh_->elements_range())
            for(SideIter side = ele.side(0); side->side_idx() < ele->n_sides(); ++side) {
                // partition
                if (init_side_ds.is_local(is)) {
                    // find (new) proc of the element of the side
                    loc_part[loc_i++] = el_ds->get_proc(
                            row_4_el[ side->element().idx() ]);
                }
                // id array
                id_4_old[is++] = side_dof( side );
            }
    }

    Partitioning::id_maps(mesh_->n_sides(), id_4_old, init_side_ds, loc_part, side_ds,
            side_id_4_loc, side_row_4_id);
    delete [] loc_part;
    delete [] id_4_old;

    // convert row_4_id arrays from separate numberings to global numbering of rows
    make_row_numberings();

}

// ========================================================================
// to finish row_4_id arrays we have to convert individual numberings of
// sides/els/edges to whole numbering of rows. To this end we count shifts
// for sides/els/edges on each proc and then we apply them on row_4_id
// arrays.
// we employ macros to avoid code redundancy
// =======================================================================
void MH_DofHandler::make_row_numberings() {
    int i, shift;
    int np = edge_ds->np();
    int edge_shift[np], el_shift[np], side_shift[np];

    int edge_n_id = mesh_->n_edges(),
            el_n_id = mesh_->n_elements(),
            side_n_id = mesh_->n_sides();

    // compute shifts on every proc
    shift = 0; // in this var. we count new starts of arrays chunks
    for (i = 0; i < np; i++) {
        side_shift[i] = shift - (side_ds->begin(i)); // subtract actual start of the chunk
        shift += side_ds->lsize(i);
        el_shift[i] = shift - (el_ds->begin(i));
        shift += el_ds->lsize(i);
        edge_shift[i] = shift - (edge_ds->begin(i));
        shift += edge_ds->lsize(i);
    }
    // apply shifts
    for (i = 0; i < side_n_id; i++) {
    	LongIdx &what = side_row_4_id[i];
        if (what >= 0)
            what += side_shift[side_ds->get_proc(what)];
    }
    for (i = 0; i < el_n_id; i++) {
    	LongIdx &what = row_4_el[i];
        if (what >= 0)
            what += el_shift[el_ds->get_proc(what)];

    }
    for (i = 0; i < edge_n_id; i++) {
    	LongIdx &what = row_4_edge[i];
        if (what >= 0)
            what += edge_shift[edge_ds->get_proc(what)];
    }
}


unsigned int MH_DofHandler::side_dof(const SideIter side) const {
    return elem_side_to_global[ side->element().idx() ][ side->side_idx() ];
}


void MH_DofHandler::set_solution( double time, double * solution) {
	OLD_ASSERT( solution != NULL, "Empty solution.\n");
    mh_solution = solution;
    time_ = time;
}

/// temporary replacement for DofHandler accessor, flux through given side
double MH_DofHandler::side_flux(const Side &side) const {
    return mh_solution[ elem_side_to_global[ side.element().idx() ][ side.side_idx() ] ];
}

/// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
double MH_DofHandler::side_scalar(const Side &side) const {
    unsigned int i_edg = side.edge_idx();
    return mh_solution[ side.mesh()->n_sides() + side.mesh()->n_elements() + i_edg ];
}


double MH_DofHandler::element_scalar( ElementAccessor<3> &ele ) const {
    return mh_solution[ mesh_->n_sides() + ele.idx() ];
}
