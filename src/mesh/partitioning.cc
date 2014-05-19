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

#include "petscao.h"


namespace IT = Input::Type;
IT::Selection Partitioning::graph_type_sel
    =IT::Selection("GraphType",
            "Different algorithms to make the sparse graph with weighted edges\n"
            "from the multidimensional mesh. Main difference is dealing with \n"
            "neighborings of elements of different dimension.")
    .add_value(any_neighboring, "any_neighboring", "Add edge for any pair of neighboring elements.")
    .add_value(any_weight_lower_dim_cuts, "any_wight_lower_dim_cuts",  "Same as before and assign higher weight to cuts of lower dimension in order to make them stick to one face.")
    .add_value(same_dimension_neighboring, "same_dimension_neghboring", "Add edge for any pair of neighboring elements of same dimension (bad for matrix multiply).")
    .close();

IT::Selection Partitioning::tool_sel
    =IT::Selection("PartTool", "Select the partitioning tool to use.")
    .add_value(PETSc, "PETSc", "Use PETSc interface to various partitioning tools.")
    .add_value(METIS, "METIS", "Use direct interface to Metis.")
    .close();

IT::Record Partitioning::input_type
    = IT::Record("Partition","Setting for various types of mesh partitioning." )
    .declare_key("tool", Partitioning::tool_sel, IT::Default("METIS"),  "Software package used for partitioning. See corresponding selection.")
    .declare_key("graph_type", Partitioning::graph_type_sel, IT::Default("any_neighboring"), "Algorithm for generating graph and its weights from a multidimensional mesh.")
    .allow_auto_conversion("graph_type") // mainly in order to allow Default value for the whole record Partition
    .close();


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
    int li, e_idx;
    unsigned int i_neigh;
    int i_s, n_s;
    F_ENTRY;

    // Add nigbouring edges only for "any_*" graph types
    bool neigh_on = ( in_.val<PartitionGraphType>("graph_type") != same_dimension_neighboring );

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
        init_el_ds_ = new Distribution(DistributionBlock(), mesh_->n_elements(), mesh_->get_comm() );  // initial distr.
        graph_ = new SparseGraphPETSC(*init_el_ds_);

        break;
    case METIS:
        init_el_ds_ = new Distribution(DistributionLocalized(), mesh_->n_elements(), mesh_->get_comm() );  // initial distr.
        graph_ = new SparseGraphMETIS(*init_el_ds_);
        break;
    }
    make_element_connection_graph();

    // compute partitioning
    loc_part_ = new int[init_el_ds_->lsize()];
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
void Partitioning::id_maps(int n_ids, int *id_4_old,
        const Distribution &old_ds, int *loc_part,
        Distribution * &new_ds, int * &id_4_loc, int * &new_4_id) {

    IS part, new_numbering;
    unsigned int size = old_ds.size(); // whole size of distr. array
    int new_counts[old_ds.np()];
    AO new_old_ao;
    int *old_4_new;
    int i_loc;
    F_ENTRY;
    // make distribution and numbering
    //DBGPRINT_INT("Local partitioning",old_ds->lsize,loc_part_);

    ISCreateGeneral(PETSC_COMM_WORLD, old_ds.lsize(), loc_part, PETSC_COPY_VALUES, &part); // global IS part.
    ISPartitioningCount(part, old_ds.np(), new_counts); // new size of each proc

    new_ds = new Distribution((unsigned int *) new_counts, PETSC_COMM_WORLD); // new distribution
    ISPartitioningToNumbering(part, &new_numbering); // new numbering

    //xprintf(Msg,"Func: %d\n",petscstack->currentsize);
    //   xprintf(Msg,"Func: %s\n",petscstack->function[petscstack->currentsize]);
    //xprintf(Msg,"Func: %s\n",petscstack->function[petscstack->currentsize-1]);

    old_4_new = (int *) xmalloc(size * sizeof(int));
    id_4_loc = (int *) xmalloc(new_ds->lsize() * sizeof(int));
    new_4_id = (int *) xmalloc((n_ids + 1) * sizeof(int));

    // create whole new->old mapping on each proc
    //DBGMSG("Creating global new->old mapping ...\n");
    AOCreateBasicIS(new_numbering, PETSC_NULL, &new_old_ao); // app ordering= new; petsc ordering = old
    for (unsigned int i = 0; i < size; i++)
        old_4_new[i] = i;
    AOApplicationToPetsc(new_old_ao, size, old_4_new);
    AODestroy(&(new_old_ao));

    // compute id_4_loc
    //DBGMSG("Creating loc.number -> id mapping ...\n");
    i_loc = 0;
    //DBGPRINT_INT("id_4_old",old_ds.lsize(),id_4_old);
    //DBGPRINT_INT("old_4_new",new_ds->lsize(),old_4_new)

    for (unsigned int i_new = new_ds->begin(); i_new < new_ds->end(); i_new++) {
        //printf("i_new: %d old: %d id: %d i_loc: %d \n",i_new,old_4_new[i_new],i_loc);
        id_4_loc[i_loc++] = id_4_old[old_4_new[i_new]];
    }
    // compute row_4_id
    //DBGMSG("Creating id -> stiffness mtx. row mapping ...\n");
    for (i_loc = 0; i_loc <= n_ids; i_loc++)
        new_4_id[i_loc] = -1; // ensure that all ids are initialized
    for (unsigned int i_new = 0; i_new < size; i_new++)
        new_4_id[id_4_old[old_4_new[i_new]]] = i_new;
    xfree(old_4_new);
}


void Partitioning::id_maps(int n_ids, int *id_4_old,  Distribution * &new_ds, int * &id_4_loc, int * &new_4_id) {
    Partitioning::id_maps(n_ids, id_4_old, *init_el_ds_, loc_part_, new_ds, id_4_loc, new_4_id);
}



vector<int> &Partitioning::seq_output_partition() {
    ASSERT(loc_part_, "Partition is not yet computed.\n");
    if (seq_part_.size() == 0) {
    //    cout << "[" << myp() << "]" << seqDBGMSG("make distr part\n");
        if (init_el_ds_->myp() == 0)
            seq_part_.resize(init_el_ds_->size());
        else
            seq_part_.resize(1);
        //for(unsigned int i = 0; i<init_el_ds_.lsize();i++) seq_part_[init_el_ds_->begin()+i]=loc
        // communicate send sizes


        MPI_Gatherv(loc_part_, init_el_ds_->lsize(), MPI_INT,
                &seq_part_[0],
                (int *)(init_el_ds_->get_lsizes_array()),
                (int *)(init_el_ds_->get_starts_array()),
                MPI_INT, 0,init_el_ds_->get_comm() );

    }
    //ASSERT_EQUAL(seq_part_.size(), init_el_ds_->size());
    return seq_part_;

}
