/*
 * partitioning.cc
 *
 *  Created on: May 3, 2013
 *      Author: jb
 */

#include "mesh/partitioning.hh"
#include "mesh/mesh.h"
#include "la/sparse_graph.hh"
#include "la/distribution.hh"

namespace IT = Input::Type;
IT::Selection Partitioning::graph_type_sel
    =IT::Selection("GraphType",
            "Different algorithms to make the sparse graph with weighted edges\n"
            "from the multidimensional mesh. Main difference is dealing with \n"
            "neighborings of elements of different dimension.")
    .add_value(any_neighboring, "any_neighboring", "Add edge for any pair of neighboring elements.")
    .add_value(any_wight_lower_dim_cuts, "any_wight_lower_dim_cuts",  "Same as before and assign higher weight to cuts of lower dimension in order to make them stick to one face.")
    .add_value(same_dimension_neghboring, "same_dimension_neghboring", "Add edge for any pair of neighboring elements of same dimension (bad for matrix multiply).")
    .close();

IT::Selection Partitioning::tool_sel
    =IT::Selection("PartTool", "Select the partitioning tool to use.")
    .add_value(PETSc, "PETSc", "Use PETSc interface to various partitioning tools.")
    .add_value(METIS, "METIS", "Use direct interface to Metis.")
    .close();

IT::Record Partitioning::input_type
    = IT::Record("Partition","Setting for various types of mesh partitioning." )
    .declare_key("tool", Partitioning::tool_sel, IT::Default("metis"),  "Software package used for partitioning. See corresponding selection.")
    .declare_key("graph_type", Partitioning::graph_type_sel, IT::Default(""), "Algorithm for generating graph and its weights from a multidimensional mesh.")
    .close();


Partitioning::Partitioning(Mesh *mesh, Input::Record in)
: mesh_(mesh), in_(in), graph_(NULL),loc_part_(NULL), init_el_ds_(NULL)
{}



Partitioning::~Partitioning() {
    if (loc_part_) delete [] loc_part_;
    loc_part_ = NULL;
    if (init_el_ds_) delete [] init_el_ds_;
    init_el_ds_ = NULL;
}

const Distribution * Partitioning::get_init_distr() const {
    ASSERT(init_el_ds_, "NULL initial distribution.");
    return init_el_ds_;
}



const int * Partitioning::get_loc_part() const {
    ASSERT(loc_part_, "NULL local partitioning.");
    return loc_part_;
}



void Partitioning::make_element_connection_graph() {

    Distribution edistr = graph_->get_distr();

    Edge *edg;
    int li, e_idx, i_neigh;
    int i_s, n_s;
    F_ENTRY;

    // Add nigbouring edges only for "any_*" graph types
    bool neigh_on = ( in_.val<PartitionGraphType>("graph_type") != same_dimension_neghboring );

    FOR_ELEMENTS( mesh_, ele) {
        //xprintf(Msg,"Element id %d , its index %d.\n",ele.id(), i_ele);

        // skip non-local elements
        if (!edistr.is_local(ele.index()))
            continue;

        // for all connected elements
        FOR_ELEMENT_SIDES( ele, si ) {
            edg = ele->side(si)->edge();

            FOR_EDGE_SIDES( edg, li ) {
                ASSERT(edg->side(li)->valid(),"NULL side of edge.");
                e_idx = ELEMENT_FULL_ITER(mesh_, edg->side(li)->element()).index();

                // for elements of connected elements, excluding element itself
                if (e_idx != ele.index()) {
                    graph_->set_edge(ele.index(), e_idx);
                }
            }
        }

        // include connections from lower dim. edge
        // to the higher dimension
        if (neigh_on) {
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb; i_neigh++) {
               n_s = ele->neigh_vb[i_neigh]->edge()->n_sides;
                for (i_s = 0; i_s < n_s; i_s++) {
                   e_idx=ELEMENT_FULL_ITER(mesh_, ele->neigh_vb[i_neigh]->edge()->side(i_s)->element()).index();
                    graph_->set_edge(ele.index(), e_idx);
                    graph_->set_edge(e_idx, ele.index());
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
        init_el_ds_ = new Distribution(Distribution::Block, mesh_->n_elements());  // initial distr.
        graph_ = new SparseGraphPETSC(*init_el_ds_);

        break;
    case METIS:
        init_el_ds_ = new Distribution(mesh_->n_elements());  // initial distr.
        graph_ = new SparseGraphMETIS(mesh_->n_elements());
        break;
    }
    make_element_connection_graph();

    // compute partitioning
    int * loc_part_ = new int[init_el_ds_->lsize()];
    graph_->partition(loc_part_);
    delete graph_;
}
