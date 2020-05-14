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
 * @file    partitioning.cc
 * @brief   
 */

#include "system/index_types.hh"
#include "mesh/partitioning.hh"
#include "la/sparse_graph.hh"
#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/neighbours.h"

#include "petscao.h"


namespace IT = Input::Type;
const IT::Selection & Partitioning::get_graph_type_sel() {
	return IT::Selection("GraphType",
            "Different algorithms to make the sparse graph with weighted edges\n"
            "from the multidimensional mesh. Main difference is dealing with \n"
            "neighboring of elements of different dimension.")
		.add_value(any_neighboring, "any_neighboring", "Add an edge for any pair of neighboring elements.")
		.add_value(any_weight_lower_dim_cuts, "any_weight_lower_dim_cuts",  "Same as before and assign higher weight to cuts of lower dimension in order to make them stick to one face.")
		.add_value(same_dimension_neighboring, "same_dimension_neighboring", "Add an edge for any pair of neighboring elements of the same dimension (bad for matrix multiply).")
		.close();
}

const IT::Selection & Partitioning::get_tool_sel() {
	return IT::Selection("PartTool", "Select the partitioning tool to use.")
		.add_value(PETSc, "PETSc", "Use PETSc interface to various partitioning tools.")
		.add_value(METIS, "METIS", "Use direct interface to Metis.")
		.close();
}

const IT::Record & Partitioning::get_input_type() {
    static IT::Record input_type = IT::Record("Partition","Setting for various types of mesh partitioning." )
		.declare_key("tool", Partitioning::get_tool_sel(), IT::Default("\"METIS\""),  "Software package used for partitioning. See corresponding selection.")
		.declare_key("graph_type", Partitioning::get_graph_type_sel(), IT::Default("\"any_neighboring\""), "Algorithm for generating graph and its weights from a multidimensional mesh.")
		.allow_auto_conversion("graph_type") // mainly in order to allow Default value for the whole record Partition
		.close();
    input_type.finish();

    return input_type;
}


Partitioning::Partitioning(Mesh *mesh, Input::Record in)
: mesh_(mesh), in_(in), graph_(NULL), loc_part_(NULL), init_el_ds_(NULL)
{
    make_partition();
}



Partitioning::~Partitioning() {
    if (loc_part_) delete [] loc_part_;
    loc_part_ = NULL;
    if (init_el_ds_) delete init_el_ds_;
    init_el_ds_ = NULL;
    if (graph_) delete graph_;
    graph_ = NULL;
}

const Distribution * Partitioning::get_init_distr() const {
	OLD_ASSERT(init_el_ds_, "NULL initial distribution.");
    return init_el_ds_;
}



const LongIdx * Partitioning::get_loc_part() const {
	OLD_ASSERT(loc_part_, "NULL local partitioning.");
    return loc_part_;
}



void Partitioning::make_element_connection_graph() {

    Distribution edistr = graph_->get_distr();

    unsigned int e_idx;
    unsigned int i_neigh;
    int i_s, n_s;

    // Add nigbouring edges only for "any_*" graph types
    bool neigh_on = ( in_.val<PartitionGraphType>("graph_type") != same_dimension_neighboring );

    for (auto ele : mesh_->elements_range()) {
        // skip non-local elements
        if ( !edistr.is_local( ele.idx() ) )
            continue;

        // for all connected elements
        for (unsigned int si=0; si<ele->n_sides(); si++) {
            Edge edg = ele.side(si)->edge();

            for (unsigned int li=0; li<edg.n_sides(); li++) {
            	ASSERT(edg.side(li)->is_valid()).error("NULL side of edge.");
                e_idx = edg.side(li)->element().idx();

                // for elements of connected elements, excluding element itself
                if ( e_idx != ele.idx() ) {
                    graph_->set_edge( ele.idx(), e_idx );
                }
            }
        }

        // include connections from lower dim. edge
        // to the higher dimension
        if (neigh_on) {
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb(); i_neigh++) {
               n_s = ele->neigh_vb[i_neigh]->edge().n_sides();
                for (i_s = 0; i_s < n_s; i_s++) {
                   e_idx = ele->neigh_vb[i_neigh]->edge().side(i_s)->element().idx();
                    graph_->set_edge( ele.idx(), e_idx );
                    graph_->set_edge( e_idx, ele.idx() );
                }
            }
        }
    }
    graph_->finalize();
}



void Partitioning::make_partition() {


    // prepare dual graph
    switch ( in_.val<PartitionTool>("tool")) {
    case PETSc:
        init_el_ds_ = new Distribution(DistributionBlock(), mesh_->n_elements(), mesh_->get_comm() );  // initial distr.
        graph_ = new SparseGraphPETSC(*init_el_ds_);

        break;
    case METIS:
        init_el_ds_ = new Distribution(DistributionLocalized(), mesh_->n_elements(), mesh_->get_comm() );  // initial distr.
        graph_ = new SparseGraphMETIS(*init_el_ds_);
        break;
    }
    int mesh_size = mesh_->n_elements();
    int num_of_procs = init_el_ds_->np();
    if (mesh_size < num_of_procs) { // check if decomposing is possible
        THROW( ExcDecomposeMesh() << EI_NElems( mesh_size ) << EI_NProcs( num_of_procs ) );
    }

    make_element_connection_graph();

    // compute partitioning
    loc_part_ = new LongIdx[init_el_ds_->lsize()];
    graph_->partition(loc_part_);
    delete graph_; graph_ = NULL;
}




/**
 * Old UGLY, PETSC dependent method for getting new numbering after partitioning.
 *
 * n_ids - given maximal ID used in id_4_old
 * id_4_old - given array of size init_el_ds_.size() - assign ID to an old index
 *
 * new_ds - new distribution of elements according to current distributed partitioning loc_part_
 * id_4_loc - IDs for local elements in new distribution, has size new_ds->lsize()
 * new_4_id - for given ID, the new index, -1 for unknown IDs
 *
 */
void Partitioning::id_maps(int n_ids, LongIdx *id_4_old,
        const Distribution &old_ds, LongIdx *loc_part,
        Distribution * &new_ds, LongIdx * &id_4_loc, LongIdx * &new_4_id) {

    IS part, new_numbering;
    unsigned int size = old_ds.size(); // whole size of distr. array
    int new_counts[old_ds.np()];
    AO new_old_ao;
    int *old_4_new;
    int i_loc;

    // make distribution and numbering
    ISCreateGeneral(PETSC_COMM_WORLD, old_ds.lsize(), loc_part, PETSC_COPY_VALUES, &part); // global IS part.
    ISPartitioningCount(part, old_ds.np(), new_counts); // new size of each proc

    new_ds = new Distribution((unsigned int *) new_counts, PETSC_COMM_WORLD); // new distribution
    ISPartitioningToNumbering(part, &new_numbering); // new numbering
    ISDestroy(&part);

    old_4_new = new int [size];
    id_4_loc = new LongIdx [ new_ds->lsize() ];
    new_4_id = new LongIdx [ n_ids + 1 ];

    // create whole new->old mapping on each proc
    AOCreateBasicIS(new_numbering, PETSC_NULL, &new_old_ao); // app ordering= new; petsc ordering = old
    ISDestroy(&new_numbering);
    for (unsigned int i = 0; i < size; i++)
        old_4_new[i] = i;
    AOApplicationToPetsc(new_old_ao, size, old_4_new);
    AODestroy(&(new_old_ao));

    // compute id_4_loc
    i_loc = 0;

    for (unsigned int i_new = new_ds->begin(); i_new < new_ds->end(); i_new++) {
        id_4_loc[i_loc++] = id_4_old[old_4_new[i_new]];
    }
    // compute row_4_id
    for (i_loc = 0; i_loc <= n_ids; i_loc++)
        new_4_id[i_loc] = -1; // ensure that all ids are initialized
    for (unsigned int i_new = 0; i_new < size; i_new++)
        new_4_id[id_4_old[old_4_new[i_new]]] = i_new;


    delete [] old_4_new;
}


void Partitioning::id_maps(int n_ids, LongIdx *id_4_old,  Distribution * &new_ds, LongIdx * &id_4_loc, LongIdx * &new_4_id) {
    Partitioning::id_maps(n_ids, id_4_old, *init_el_ds_, loc_part_, new_ds, id_4_loc, new_4_id);
}



shared_ptr< vector<int> > Partitioning::subdomain_id_field_data() {
    ASSERT(loc_part_).error("Partition is not yet computed.\n");
    if (!seq_part_) {
    	unsigned int seq_size=mesh_->get_el_ds()->lsize();
    	seq_part_ = make_shared< vector<int> >(seq_size, mesh_->get_el_ds()->myp());
    }

    return seq_part_;

}
