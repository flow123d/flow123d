/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Setup and solve linear system of mixed-hybrid discretization of the linear
 * porous media flow with possible preferential flow in fractures and chanels.
 *
 *
 *
 *
 */

#include "petscmat.h"
#include "petscviewer.h"
#include "petscao.h"
#include "petscerror.h"

#include "system.hh"
#include "math_fce.h"
#include "mesh.h"
#include "par_distribution.hh"
#include "darcy_flow_mh.hh"
#include "la_linsys.hh"
#include "solve.h"
#include "la_schur.hh"
#include "sparse_graph.hh"

//static void compute_nonzeros( DarcyFlowMH *w);
static void make_schur0(DarcyFlowMH *w);
static void make_schur1(DarcyFlowMH *w);
static void make_schur2(DarcyFlowMH *w);

static void destroy_water_linsys(DarcyFlowMH *w);

//#if USE_PARALLEL
//#include "parmetis.h"


//void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *);
//#endif

// PRIVATE HEADERS
static void prepare_parallel(DarcyFlowMH *w);

//=============================================================================
// CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
// - do it in parallel:
//   - initial distribution of elements, edges
//
/*! @brief CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
 *
 * Parameters {Solver,NSchurs} number of performed Schur
 * complements (0,1,2) for water flow MH-system
 *
 */
//=============================================================================
void create_water_linsys(Mesh *mesh, DarcyFlowMH **XWSys) {
    DarcyFlowMH *WSys;
    int ierr;

    ASSERT(!( mesh == NULL ),"NULL mesh pointer.\n");

    // create and init new DarcyFlowMH
    if (*XWSys) {
        //TODO: jak se tohle pouziva? kdyz to odkomentuju, tak to vybouchne v PETSc
        // je nutne odalokovat i vsechny vnorene struktury, jinak to je mem leak jako prase...

        //destroy_water_linsys(*XWSys);
        xfree(*XWSys);
    }
    *XWSys = (DarcyFlowMH*) xmalloc(sizeof(DarcyFlowMH));

    WSys = *XWSys;
    WSys->mesh = mesh;
    WSys->size = mesh->n_elements() + mesh->n_sides + mesh->n_edges;
    WSys->n_schur_compls = OptGetInt("Solver", "NSchurs", "2");
    if ((unsigned int) WSys->n_schur_compls > 2) {
        xprintf(Warn,"Invalid number of Schur Complements. Using 2.");
        WSys->n_schur_compls = 2;
    }

    WSys->solver = (Solver *) xmalloc(sizeof(Solver));
    solver_init(WSys->solver);

    WSys->solution = NULL;
    WSys->schur0 = NULL;
    WSys->schur1 = NULL;
    WSys->schur2 = NULL;

    // TODO PARALLEL
    // init paralel structure

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &(WSys->myp));
    ierr += MPI_Comm_size(PETSC_COMM_WORLD, &(WSys->np));
    if (ierr)
        xprintf(Err, "Some error in MPI.\n");
    MPI_Barrier(PETSC_COMM_WORLD);
    prepare_parallel(WSys);

    WSys->side_ds->view();
    WSys->el_ds->view();
    WSys->edge_ds->view();
    MPI_Barrier(PETSC_COMM_WORLD);

    make_schur0(WSys);
}

// ================================================
// PARALLLEL PART
//

/**
 * Make connectivity graph of the second Schur complement and compute optimal partitioning.
 * This verison assign 1D and 2D edges to one processor, represent them as
 * a weighted vertex, and 2D-3D neghbourings as weighted edges. 3D edges are
 * then distributed over all procesors.
 */
void make_edge_conection_graph(Mesh *mesh, SparseGraph * &graph) {

    Distribution edistr=graph->get_distr();
    Edge *edg;
    Element *ele;
    int li, si, eid, i_neigh, i_edg;
    int e_weight;

    int edge_dim_weights[3] = { 100, 10, 1 };
    F_ENTRY;

    i_edg = 0;
    FOR_EDGES(edg) {
        ASSERT( edg->id == i_edg, "Edge id %d doesn't match its index %d.\n",edg->id, i_edg);

        // skip non-local edges
        if (!edistr.is_local(edg->id)) {
            i_edg++;
            continue;
        }

        e_weight = edge_dim_weights[edg->side[0]->element->dim - 1];
        // for all connected elements
        FOR_EDGE_SIDES( edg, li ) {
            ASSERT(NONULL(edg->side[li]),"NULL side of edge.");
            ele = edg->side[li]->element;
            ASSERT(NONULL(ele),"NULL element of side.");

            // for sides of connected element, excluding edge itself
            FOR_ELEMENT_SIDES( ele, si ) {
                eid = ele->side[si]->edge->id;
                if (eid != edg->id)
                    graph->set_edge(edg->id, eid, e_weight);
            }

            // include connections from lower dim. edge
            // to the higher dimension
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb; i_neigh++) {
                eid = ele->neigh_vb[i_neigh]->edge->id;
                graph->set_edge(edg->id, eid, e_weight);
                graph->set_edge(eid, edg->id, e_weight);
            }
        }
        i_edg++;
    }

    graph->finalize();
}

/**
 * Make connectivity graph of elements of mesh - dual graph: elements vertices of graph.
 * Jakub S.
 */
void make_element_connection_graph(Mesh *mesh, SparseGraph * &graph,bool neigh_on) {

    Distribution edistr=graph->get_distr();

    Edge *edg;
    int li, si, e_idx, i_neigh;
    int i_s, n_s;
    F_ENTRY;

    FOR_ELEMENTS(ele) {
        //xprintf(Msg,"Element id %d , its index %d.\n",ele.id(), i_ele);

        // skip non-local elements
        if (!edistr.is_local(ele.index()))
            continue;

        // for all connected elements
        FOR_ELEMENT_SIDES( ele, si ) {
            edg = ele->side[si]->edge;

            FOR_EDGE_SIDES( edg, li ) {
                ASSERT(NONULL(edg->side[li]),"NULL side of edge.");
                e_idx = ELEMENT_FULL_ITER(edg->side[li]->element).index();

                // for elements of connected elements, excluding element itself
                if (e_idx != ele.index()) {
                    graph->set_edge(ele.index(), e_idx);
                }
            }
        }

        // TODO: Suitable way to represent connections between dimensions in graph.
        //
        // include connections from lower dim. edge
        // to the higher dimension
        if ( neigh_on ) {
            for(i_neigh=0; i_neigh < ele->n_neighs_vb; i_neigh++) {
               n_s = ele->neigh_vb[i_neigh]->edge->n_sides;
               for(i_s=0; i_s < n_s; i_s++) {
                   e_idx=ELEMENT_FULL_ITER(ele->neigh_vb[i_neigh]->edge->side[i_s]->element).index();
                   graph->set_edge(ele.index(),e_idx);
                   graph->set_edge(e_idx,ele.index());
               }
            }
        }
    }
    graph->finalize();
}

// ==========================================================================
// from old->id mapping and IS parttioning create:
// - new distribution, new numbering
// - loc->id (array od ids of local elements)
// - id->new (new index for each id)
// ==========================================================================
void id_maps(int n_ids, int *id_4_old, const Distribution &old_ds,
        int *loc_part, Distribution * &new_ds, int * &id_4_loc, int * &new_4_id) {
    IS part, new_numbering;
    unsigned int size = old_ds.size(); // whole size of distr. array
    int new_counts[old_ds.np()];
    AO new_old_ao;
    int *old_4_new;
    int i, i_loc, i_new, i_old;
    F_ENTRY;
    // make distribution and numbering
    //DBGPRINT_INT("Local partitioning",old_ds->lsize,loc_part);

    ISCreateGeneral(PETSC_COMM_WORLD, old_ds.lsize(), loc_part, &part); // global IS part.
    ISPartitioningCount(part, old_ds.np(), new_counts); // new size of each proc

    new_ds = new Distribution((unsigned int *) new_counts); // new distribution
    ISPartitioningToNumbering(part, &new_numbering); // new numbering

    //xprintf(Msg,"Func: %d\n",petscstack->currentsize);
    //   xprintf(Msg,"Func: %s\n",petscstack->function[petscstack->currentsize]);
    //xprintf(Msg,"Func: %s\n",petscstack->function[petscstack->currentsize-1]);

    old_4_new = (int *) xmalloc(size * sizeof(int));
    id_4_loc = (int *) xmalloc(new_ds->lsize() * sizeof(int));
    new_4_id = (int *) xmalloc((n_ids + 1) * sizeof(int));

    // create whole new->old mapping on each proc
    DBGMSG("Creating global new->old mapping ...\n");
    AOCreateBasicIS(new_numbering, PETSC_NULL, &new_old_ao); // app ordering= new; petsc ordering = old
    for (i = 0; i < size; i++)
        old_4_new[i] = i;
    AOApplicationToPetsc(new_old_ao, size, old_4_new);
    AODestroy(new_old_ao);

    // compute id_4_loc
    DBGMSG("Creating loc.number -> id mapping ...\n");
    i_loc = 0;
    //DBGPRINT_INT("id_4_old",old_ds.lsize(),id_4_old);
    //DBGPRINT_INT("old_4_new",new_ds->lsize(),old_4_new)

    for (i_new = new_ds->begin(); i_new < new_ds->end(); i_new++) {
        //printf("i_new: %d old: %d id: %d i_loc: %d \n",i_new,old_4_new[i_new],i_loc);
        id_4_loc[i_loc++] = id_4_old[old_4_new[i_new]];
    }
    // compute row_4_id
    DBGMSG("Creating id -> stiffness mtx. row mapping ...\n");
    for (i_loc = 0; i_loc <= n_ids; i_loc++)
        new_4_id[i_loc] = -1; // ensure that all ids are initialized
    for (i_new = 0; i_new < size; i_new++)
        new_4_id[id_4_old[old_4_new[i_new]]] = i_new;
    xfree(old_4_new);
}

// ========================================================================
// to finish row_4_id arrays we have to convert individual numberings of
// sides/els/edges to whole numbering of rows. To this end we count shifts
// for sides/els/edges on each proc and then we apply them on row_4_id
// arrays.
// we employ macros to avoid code redundancy
// =======================================================================
void make_row_numberings(DarcyFlowMH *w) {
    int i, shift;
    int np = w->edge_ds->np();
    int edge_shift[np], el_shift[np], side_shift[np];
    unsigned int rows_starts[np];
    Mesh *mesh = w->mesh;
    int edge_n_id = mesh->max_edg_id + 1, el_n_id = mesh->element.size(),
            side_n_id = mesh->max_side_id + 1;

    // compute shifts on every proc
    shift = 0; // in this var. we count new starts of arrays chunks
    for (i = 0; i < np; i++) {
        side_shift[i] = shift - (w->side_ds->begin(i)); // subtract actual start of the chunk
        shift += w->side_ds->lsize(i);
        el_shift[i] = shift - (w->el_ds->begin(i));
        shift += w->el_ds->lsize(i);
        edge_shift[i] = shift - (w->edge_ds->begin(i));
        shift += w->edge_ds->lsize(i);
        rows_starts[i] = shift;
    }
    //DBGPRINT_INT("side_shift",np,side_shift);
    //DBGPRINT_INT("el_shift",np,el_shift);
    //DBGPRINT_INT("edge_shift",np,edge_shift);
    // apply shifts
    for (i = 0; i < side_n_id; i++) {
        int &what = w->side_row_4_id[i];
        if (what >= 0)
            what += side_shift[w->side_ds->get_proc(what)];
    }
    for (i = 0; i < el_n_id; i++) {
        int &what = w->row_4_el[i];
        if (what >= 0)
            what += el_shift[w->el_ds->get_proc(what)];

    }
    for (i = 0; i < edge_n_id; i++) {
        int &what = w->edge_row_4_id[i];
        if (what >= 0)
            what += edge_shift[w->edge_ds->get_proc(what)];
    }
    // make distribution of rows
    for (i = np - 1; i > 0; i--)
        rows_starts[i] -= rows_starts[i - 1];
    w->rows_ds = new Distribution(rows_starts);
}

// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void prepare_parallel(DarcyFlowMH *w) {

    int *loc_part; // optimal (edge,el) partitioning (local chunk)
    int *id_4_old; // map from old idx to ids (edge,el)
    // auxiliary
    Mesh *mesh = w->mesh;
    Edge *edg;
    Element *el;
    Side *side;
    int i, loc_i, x_proc;

    int n_edg, n_e, n_sides, ndof, idof;
    int n_s, i_s, ig4s;
    int i_neigh;
    int sid, edgid, e_idx;
    int lmap_aux;
    int *map_aux;
    int myid;
    int ind_row;
    int i_loc, el_row, side_row, edge_row, nsides;

    int ndof_loc;
    PetscErrorCode err;
    Mat sub_matrix;
    F_ENTRY;
    MPI_Barrier(PETSC_COMM_WORLD);

    if (w->solver->type == PETSC_MATIS_SOLVER) {
        xprintf(Msg,"Compute optimal partitioning of elements.\n");

        // prepare dual graph
        Distribution init_el_ds(Distribution::Block, mesh->n_elements());  // initial distr.
        SparseGraph *element_graph= new SparseGraphPETSC(init_el_ds);
        int *loc_part = new int[(init_el_ds.size() + 1)];                                     // partitionig in initial distribution

        make_element_connection_graph(mesh, element_graph);
        WARN_ASSERT(element_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

        element_graph->partition(loc_part);

        MPI_Bcast(loc_part, init_el_ds.size(), MPI_INT, 0, MPI_COMM_WORLD);
        DBGPRINT_INT("loc_part",init_el_ds.size(),loc_part);

        // prepare parallel distribution of dofs linked to elements
        id_4_old = (int *) xmalloc(mesh->n_elements() * sizeof(int));
        i = 0;
        FOR_ELEMENTS(el)
            id_4_old[i++] = el.index();
        id_maps(mesh->element.size(), id_4_old, init_el_ds, loc_part, w->el_ds,
                w->el_4_loc, w->row_4_el);
        free(loc_part);
        free(id_4_old);

        DBGMSG("Compute appropriate edge partitioning ...\n");
        //optimal element part; loc. els. id-> new el. numbering
        Distribution init_edge_ds(Distribution::Localized, mesh->n_edges);
        // partitioning of edges, edge belongs to the proc of his first element
        // this is not optimal but simple
        loc_part = (int *) xmalloc((init_edge_ds.lsize() + 1) * sizeof(int));
        id_4_old = (int *) xmalloc(mesh->n_edges * sizeof(int));
        {
            int iedg = 0;
            loc_i = 0;
            FOR_EDGES( edg ) {
                // partition
                edgid = edg->id;
                e_idx = mesh->element.index(edg->side[0]->element);
                //xprintf(Msg,"Index of edge: %d first element: %d \n",edgid,e_idx);
                if (init_edge_ds.is_local(iedg)) {
                    // find (new) proc of the first element of the edge
                    loc_part[loc_i++] = w->el_ds->get_proc(w->row_4_el[e_idx]);
                }
                // id array
                id_4_old[iedg] = edgid;
                iedg = iedg + 1;
            }
        }
        //    // make trivial part
        //    for(loc_i=0;loc_i<init_el_ds->lsize;loc_i++) loc_part[loc_i]=init_el_ds->myp;
        //DBGPRINT_INT("loc_part",init_edge_ds.lsize(),loc_part);

        id_maps(mesh->max_edg_id, id_4_old, init_edge_ds, loc_part, w->edge_ds,
                w->edge_id_4_loc, w->edge_row_4_id);
        free(loc_part);
        free(id_4_old);

    } else {
        xprintf(Msg,"Compute optimal partitioning of edges.\n");

        SparseGraph *edge_graph = new SparseGraphMETIS(mesh->n_edges);                     // graph for partitioning
        Distribution init_edge_ds = edge_graph->get_distr();  // initial distr.
        int *loc_part = new int[init_edge_ds.lsize()];                                     // partitionig in initial distribution

        make_edge_conection_graph(mesh, edge_graph);
        WARN_ASSERT(edge_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

        edge_graph->partition(loc_part);

        delete edge_graph;


        // debugging output
/*
        if (init_edge_ds.myp() == 0) {
            Edge *edg;
            int i_edg = 0;
            int stat[3][init_edge_ds.np()];
            for (int ip = 0; ip < init_edge_ds.np(); ip++) {
                stat[0][ip] = stat[1][ip] = stat[2][ip] = 0;
            }
            for(i_edg=0;i_edg < ;i_edg++) {
                DBGMSG("edg: %d %d %d\n",
                       i_edg,edg->side[0]->element->dim-1,loc_part[i_edg]);
                int dim=edg->side[0]->element->dim - 1;
                int part=loc_part[i_edg];
                (stat[dim][part])++;
                i_edg++;
            }
            for (int ip = 0; ip < init_edge_ds.np(); ip++) {
                DBGMSG("1D: %10d 2d: %10d 3d: %10d\n",
                        stat[0][ip],stat[1][ip],stat[2][ip]);
            }
        }
*/
        id_4_old = (int *) xmalloc(mesh->n_edges * sizeof(int));
        i = 0;
        FOR_EDGES(edg)
            id_4_old[i++] = edg->id;
        id_maps(mesh->max_edg_id, id_4_old, init_edge_ds, (int *) loc_part,
                w->edge_ds, w->edge_id_4_loc, w->edge_row_4_id);


        delete loc_part;
        xfree(id_4_old);

        DBGMSG("Compute appropriate element partitioning ...\n");
        //optimal element part; loc. els. id-> new el. numbering
        Distribution init_el_ds(Distribution::Block, mesh->n_elements());
        // partitioning of elements, element belongs to the proc of his first edge
        // this is not optimal but simple
        loc_part = (int *) xmalloc(init_el_ds.lsize() * sizeof(int));
        id_4_old = (int *) xmalloc(mesh->n_elements() * sizeof(int));
        {
            int iel = 0, i_edg;
            loc_i = 0;
            FOR_ELEMENTS( el ) {
                // partition
                if (init_el_ds.is_local(iel)) {
                    // find (new) proc of the first edge of element
                    //DBGMSG("%d %d %d %d\n",iel,loc_i,el->side[0]->edge->id,w->edge_row_4_id[el->side[0]->edge->id]);
                    loc_part[loc_i++] = w->edge_ds->get_proc(
                            el->side[0]->edge->id);
                }
                // id array
                id_4_old[iel++] = mesh->element.index(el);
            }
        }
        //    // make trivial part
        //    for(loc_i=0;loc_i<init_el_ds->lsize;loc_i++) loc_part[loc_i]=init_el_ds->myp;
        id_maps(mesh->element.size(), id_4_old, init_el_ds, loc_part, w->el_ds,
                w->el_4_loc, w->row_4_el);
        xfree(loc_part);
        xfree(id_4_old);
    }

    DBGMSG("Compute side partitioning ...\n");
    //optimal side part; loc. sides; id-> new side numbering
    Distribution init_side_ds(Distribution::Block, mesh->n_sides);
    // partitioning of sides follows elements
    loc_part = (int *) xmalloc(init_side_ds.lsize() * sizeof(int) + 1);
    id_4_old = (int *) xmalloc(mesh->n_sides * sizeof(int));
    {
        int is = 0, iel;
        loc_i = 0;
        FOR_SIDES( side ) {
            // partition
            if (init_side_ds.is_local(is)) {
                // find (new) proc of the element of the side
                loc_part[loc_i++] = w->el_ds->get_proc(
                        w->row_4_el[mesh->element.index(side->element)]);
            }
            // id array
            id_4_old[is++] = side->id;
        }
    }
    // make trivial part
    //for(loc_i=0;loc_i<init_side_ds->lsize;loc_i++) loc_part[loc_i]=init_side_ds->myp;

    id_maps(mesh->max_side_id, id_4_old, init_side_ds, loc_part, w->side_ds,
            w->side_id_4_loc, w->side_row_4_id);
    xfree(loc_part);
    xfree(id_4_old);

    /*
     DBGPRINT_INT("edge_id_4_loc",w->edge_ds->lsize,w->edge_id_4_loc);
     DBGPRINT_INT("el_4_loc",w->el_ds->lsize,w->el_4_loc);
     DBGPRINT_INT("side_id_4_loc",w->side_ds->lsize,w->side_id_4_loc);
     DBGPRINT_INT("edge_row_4_id",mesh->n_edges,w->edge_row_4_id);
     DBGPRINT_INT("el_row_4_id",mesh->max_elm_id+1,w->el_row_4_id);
     DBGPRINT_INT("side_row_4_id",mesh->max_side_id+1,w->side_row_4_id);
     */
    // convert row_4_id arrays from separate numberings to global numbering of rows
    //MPI_Barrier(PETSC_COMM_WORLD);
    //DBGMSG("Finishing row_4_id\n");
    //MPI_Barrier(PETSC_COMM_WORLD);
    make_row_numberings(w);
    //DBGPRINT_INT("edge_row_4_id",mesh->n_edges,w->edge_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh->max_elm_id+1,w->el_row_4_id);
    //DBGPRINT_INT("side_row_4_id",mesh->max_side_id+1,w->side_row_4_id);

    w->lsize = w->side_ds->lsize() + w->el_ds->lsize() + w->edge_ds->lsize();

    // make old_4_new
    w->old_4_new = (int *) malloc((mesh->n_edges + mesh->n_sides
            + mesh->n_elements()) * sizeof(int));
    i = 0;
    FOR_SIDES( side )
        w->old_4_new[w->side_row_4_id[side->id]] = i++;
    FOR_ELEMENTS( el )
        w->old_4_new[w->row_4_el[el.index()]] = i++;
    FOR_EDGES(edg)
        w->old_4_new[w->edge_row_4_id[edg->id]] = i++;

    // prepare global_row_4_sub_row
    if (w->solver->type == PETSC_MATIS_SOLVER) {
        xprintf(Msg,"Compute mapping of local subdomain rows to global rows.\n");

        // prepare arrays of velocities, pressures and Lagrange multipliers
        n_edg = mesh->n_edges;
        n_e = mesh->n_elements();
        n_sides = mesh->n_sides;

        ndof = n_edg + n_e + n_sides;
        xprintf(Msg,"n_edg = %d n_e = %d n_sides = %d ndof = %d \n ",n_edg,n_e,n_sides,ndof);
        // TODO: use quick sort and short arrays
        // initialize array
        lmap_aux = ndof;
        map_aux = (int *) xmalloc(lmap_aux * sizeof(int) + 1);
        for (idof = 0; idof < ndof; idof++) {
            map_aux[idof] = 0;
        }

        // ordering of dofs
        // for each subdomain:
        // | velocities (at sides) | pressures (at elements) | L. mult. (at edges) |
        //
        //DBGPRINT_INT("el_4_loc",w->el_ds->lsize(),w->el_4_loc);

        // processor ID
        myid = w->el_ds->myp();

        for (i_loc = 0; i_loc < w->el_ds->lsize(); i_loc++) {
            el = mesh->element(w->el_4_loc[i_loc]);
            el_row = w->row_4_el[w->el_4_loc[i_loc]];

            map_aux[el_row] = map_aux[el_row] + 1;

            nsides = el->n_sides;
            for (i = 0; i < nsides; i++) {
                side_row = w->side_row_4_id[el->side[i]->id];
                edge_row = w->edge_row_4_id[el->side[i]->edge->id];

                map_aux[side_row] = map_aux[side_row] + 1;
                map_aux[edge_row] = map_aux[edge_row] + 1;
                xprintf(Msg,"el_row %d side_row = %d edge_row = %d \n ",el_row,side_row,edge_row);
            }

            for (i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
                edgid = el->neigh_vb[i_neigh]->edge->id;
                // mark this edge at map_aux
                edge_row = w->edge_row_4_id[edgid];
                map_aux[edge_row] = map_aux[edge_row] + 1;
                xprintf(Msg,"el_row %d edge_row = %d \n ",el_row,edge_row);
            }
        }
        //DBGPRINT_INT("map_aux",lmap_aux,map_aux);

        // count nonzeros in map_aux
        ndof_loc = 0;
        for (i = 0; i < lmap_aux; i++) {
            if (map_aux[i] > 0) {
                ndof_loc = ndof_loc + 1;
            }
        }
        xprintf(Msg,"ndof_loc = %d \n",ndof_loc);

        // initialize mapping arrays in MATIS matrix
        w->ndof_loc = ndof_loc;
        w->global_row_4_sub_row = (int *) xmalloc((ndof_loc + 1) * sizeof(int));

        ig4s = 0;
        for (i = 0; i < lmap_aux; i++) {
            if (map_aux[i] > 0) {
                w->global_row_4_sub_row[ig4s] = i;
                ig4s = ig4s + 1;
            }
        }
        // check that the array was filled
        if (ig4s != ndof_loc) {
            xprintf(PrgErr,"Data length mismatch! %d not like %d \n",ndof_loc,ig4s);
        }

        free(map_aux);

        DBGPRINT_INT("global_row_4_sub_row",w->ndof_loc,w->global_row_4_sub_row);

    }
}

void mat_count_off_proc_values(Mat m, Vec v) {
    int n,first,last;
    const PetscInt *cols;
    Distribution distr(v);

    int n_off=0;
    int n_on=0;
    int n_off_rows=0;
    MatGetOwnershipRange(m,&first,&last);
    for(int row=first;row<last;row++) {
        MatGetRow(m,row,&n,&cols,PETSC_NULL);
        bool exists_off=false;
        for(int i=0;i<n;i++)
            if (distr.get_proc(cols[i]) != distr.myp() ) n_off++,exists_off=true;
            else n_on++;
        if (exists_off) n_off_rows++;
        MatRestoreRow(m,row,&n,&cols,PETSC_NULL);
    }
    printf("[%d] rows: %d off_rows: %d on: %d off: %d\n",distr.myp(),last-first,n_off_rows,n_on,n_off);
}

// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================


// ******************************************
// ABSTRACT ASSEMBLY OF MH matrix
// TODO: matice by se mela sestavovat zvlast pro kazdou dimenzi (objem, pukliny, pruseciky puklin)
//       konekce by se mely sestavovat cyklem pres konekce, konekce by mely byt paralelizovany podle
//       distribuce elementu nizssi dimenze
//       k tomuto je treba nejdriv spojit s JK verzi, aby se vedelo co se deje v transportu a
//       predelat mesh a neigbouring
// *****************************************

void mh_abstract_assembly(DarcyFlowMH *w) {
    LinSys *ls = w->schur0;
    Mesh* mesh = w->mesh;
    ElementFullIter ele = ELEMENT_FULL_ITER(NULL);
    struct Edge *edg;

    int el_row, side_row, edge_row;
    int tmp_rows[100];
    int i, i_loc, nsides, li, si;
    int side_rows[4], edge_rows[4]; // rows for sides and edges of one element
    double f_val;
    double zeros[1000]; // to make space for second schur complement, max. 10 neigbour edges of one el.
    double minus_ones[4] = { -1.0, -1.0, -1.0, -1.0 };
    F_ENTRY;

    //DBGPRINT_INT("side_row_4_id",mesh->max_side_id+1,w->side_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh->max_elm_id+1,w->el_row_4_id);
    //DBGPRINT_INT("edge_row_4_id",mesh->max_edg_id+1,w->edge_row_4_id);
    //DBGPRINT_INT("el_id_4_loc",w->el_ds->lsize(),w->el_id_4_loc);

    SET_ARRAY_ZERO(zeros,1000);
    for (i_loc = 0; i_loc < w->el_ds->lsize(); i_loc++) {
        ele = mesh->element(w->el_4_loc[i_loc]);
        el_row = w->row_4_el[w->el_4_loc[i_loc]];
        nsides = ele->n_sides;
        for (i = 0; i < nsides; i++) {
            side_row = side_rows[i] = w->side_row_4_id[ele->side[i]->id];
            edge_row = edge_rows[i] = w->edge_row_4_id[ele->side[i]->edge->id];
            // set block C and C': side-edge, edge-side
            ls->mat_set_value(side_row, edge_row, ele->side[i]->c_val);
            ls->mat_set_value(edge_row, side_row, ele->side[i]->c_val);
        }
        // set block A: side-side on one element - block diagonal matrix
        ls->mat_set_values(nsides, side_rows, nsides, side_rows, ele->loc);
        // set block B, B': element-side, side-element
        ls->mat_set_values(1, &el_row, nsides, side_rows, minus_ones);
        ls->mat_set_values(nsides, side_rows, 1, &el_row, minus_ones);
        // set RHS for sides - dirichlet BC; RHS for elements - neuman BC
        ls->rhs_set_values(nsides, side_rows, ele->rhs);
        ls->rhs_set_value(el_row, ele->rhs_b);

        // D block: non-compatible conections and diagonal: element-element
        for (i = 0; i < ele->d_row_count; i++)
            tmp_rows[i] = w->row_4_el[ele->d_el[i]];
        ls->mat_set_values(1, &el_row, ele->d_row_count, tmp_rows, ele->d_val);
        // E',E block: compatible connections: element-edge
        for (i = 0; i < ele->e_row_count; i++)
            tmp_rows[i] = w->edge_row_4_id[ele->e_edge_id[i]];
        ls->mat_set_values(1, &el_row, ele->e_row_count, tmp_rows, ele->e_val);
        ls->mat_set_values(ele->e_row_count, tmp_rows, 1, &el_row, ele->e_val);

        // add virtual values for schur complement allocation
        switch (w->n_schur_compls) {
        case 2:
            if (ele->d_row_count > 1) {
                xprintf(Warn,"Can not use second Schur complement for problem with non-compatible connections.\n");
                w->n_schur_compls = 1;
            }
            // for 2. Schur: N dim edge is conected with N dim element =>
            // there are nz between N dim edge and N-1 dim edges of the element
            ASSERT(ele->e_row_count*nsides<1000,"Too many values in E block.");
            ls->mat_set_values(nsides, edge_rows, ele->e_row_count, tmp_rows,
                    zeros);
            ls->mat_set_values(ele->e_row_count, tmp_rows, nsides, edge_rows,
                    zeros);
            ASSERT(ele->e_row_count*ele->e_row_count<1000,"Too many values in E block.");
            ls->mat_set_values(ele->e_row_count, tmp_rows, ele->e_row_count,
                    tmp_rows, zeros);
        case 1: // included also for case 2
            // -(C')*(A-)*B block and its transpose conect edge with its elements
            ls->mat_set_values(1, &el_row, nsides, edge_rows, zeros);
            ls->mat_set_values(nsides, edge_rows, 1, &el_row, zeros);
            // -(C')*(A-)*C block conect all edges of every element
            ls->mat_set_values(nsides, edge_rows, nsides, edge_rows, zeros);
        }
    }
    //if (! mtx->ins_mod == ALLOCATE ) {
    //    MatAssemblyBegin(mtx->A,MAT_FINAL_ASSEMBLY);
    //    MatAssemblyEnd(mtx->A,MAT_FINAL_ASSEMBLY);
    // }
    // set block F - diagonal: edge-edge from Newton BC
    for (i_loc = 0; i_loc < w->edge_ds->lsize(); i_loc++) {
        edg = mesh->edge_hash[w->edge_id_4_loc[i_loc]];
        edge_row = w->edge_row_4_id[edg->id];
        //xprintf(Msg,"F: %d %f\n",w->old_4_new[edge_row],edg->f_val);
        ls->mat_set_value(edge_row, edge_row, edg->f_val);
        ls->rhs_set_value(edge_row, edg->f_rhs);
    }
}

//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void solve_water_linsys(DarcyFlowMH *w) {

    Timing *solver_time = timing_create("SOLVING MH SYSTEM", PETSC_COMM_WORLD);
    F_ENTRY;

    printf("[%d] before solve\n",w->el_ds->myp());
    switch (w->n_schur_compls) {
    case 0: /* none */
        solve_system(w->solver, w->schur0);
        break;
    case 1: /* first schur complement of A block */
        make_schur1(w);
        solve_system(w->solver, w->schur1->get_system());
        w->schur1->resolve();
        break;
    case 2: /* second shur complement of the max. dimension elements in B block */
        make_schur1(w);
        printf("[%d] a s1\n",w->el_ds->myp());
        make_schur2(w);
        printf("[%d] as2\n",w->el_ds->myp());

        mat_count_off_proc_values(w->schur2->get_system()->get_matrix(),w->schur2->get_system()->get_solution());
        solve_system(w->solver, w->schur2->get_system());
        printf("[%d] asolve\n",w->el_ds->myp());
        w->schur2->resolve();
        printf("[%d] a rs1\n",w->el_ds->myp());
        w->schur1->resolve();
        printf("[%d] a rs2\n",w->el_ds->myp());


/*
        // experiment

        Vec tmp1,tmp2;
        VecDuplicate(w->schur2->get_system()->get_solution(),&tmp1);
        VecDuplicate(w->schur2->get_system()->get_solution(),&tmp2);
        VecCopy(w->schur2->get_system()->get_solution(),tmp1);
        for(int i=1;i<100;i++) {
            printf("it: %d\n",i);
            MatMultAdd(w->schur2->get_system()->get_matrix(),tmp1,tmp1,tmp2);
            VecSwap(tmp1,tmp2);
        }
        VecDestroy(tmp1);
        VecDestroy(tmp2);
*/
        break;
    }
    // TODO PARALLEL
    xprintf(MsgVerb,"Scattering solution vector to all processors ...\n");
    // scatter solution to all procs
    if (w->solution == NULL) {
        IS is_par, is_loc;
        int i, si, *loc_idx;
        Element *ele;
        Edge *edg;
        Mesh *mesh = w->mesh;

        // create local solution vector
        w->solution = (double *) xmalloc(w->size * sizeof(double));
        VecCreateSeqWithArray(PETSC_COMM_SELF, w->size, w->solution,
                &(w->sol_vec));

        // create seq. IS to scatter par solutin to seq. vec. in original order
        // use essentialy row_4_id arrays
        loc_idx = (int *) xmalloc(w->size * sizeof(int));
        i = 0;
        FOR_ELEMENTS(ele) {
            FOR_ELEMENT_SIDES(ele,si) {
                loc_idx[i++] = w->side_row_4_id[ele->side[si]->id];
            }
        }
        FOR_ELEMENTS(ele) {
            loc_idx[i++] = w->row_4_el[ele.index()];
        }
        FOR_EDGES(edg) {
            loc_idx[i++] = w->edge_row_4_id[edg->id];
        }
        ASSERT( i==w->size,"Size of array does not match number of fills.\n");
        //DBGPRINT_INT("loc_idx",w->size,loc_idx);
        ISCreateGeneral(PETSC_COMM_SELF, w->size, loc_idx, &(is_loc));
        xfree(loc_idx);
        VecScatterCreate(w->schur0->get_solution(), is_loc, w->sol_vec,
                PETSC_NULL, &(w->par_to_all));
        ISDestroy(is_loc);
    }
    VecScatterBegin(w->par_to_all, w->schur0->get_solution(), w->sol_vec,
            INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(w->par_to_all, w->schur0->get_solution(), w->sol_vec,
            INSERT_VALUES, SCATTER_FORWARD);

    timing_destroy(solver_time);
}

//=============================================================================
// COUMPUTE NONZEROES IN THE WATER MH MATRIX
//=============================================================================
/*
 * void compute_nonzeros( TWaterLinSys *w_ls) {

 ElementIter ele;
 struct Edge *edg;
 int row,side,i;
 Mesh* mesh = w_ls->water_eq;

 xprintf( Msg, "Computing nonzero values ...\n");
 w_ls->nonzeros=(int *)xmalloc(sizeof(int)*w_ls->size);
 row=0;
 FOR_ELEMENTS( ele ) { // count A, B', C'
 // n_sides in A, 1 in B', 1 in C'
 for(side=0;side<ele->n_sides;side++)
 w_ls->nonzeros[row++]=ele->n_sides+2;
 }
 FOR_ELEMENTS( ele ) { // count B, D, E'
 // n_sides in B, # D, # E
 w_ls->nonzeros[row++]=ele->n_sides+ele->d_row_count+ele->e_row_count;
 }
 FOR_EDGES ( edg ) { // count C, E', F
 // C - n_sides(1 on BC, 2 inside), E' - possible ngh. F-diagonal
 w_ls->nonzeros[row++]=(edg->n_sides)+((edg->neigh_vb!=NULL)?1:0)+1;
 }

 // count additional space for the valueas of schur complement
 if (w_ls->n_schur_compls > 0) {
 row=w_ls->sizeA;
 // -B'*A-*B block is diagonal and already counted
 // -B'*A-*C block conect element with its edges
 FOR_ELEMENTS( ele ) {
 w_ls->nonzeros[row++]+=ele->n_sides;
 }

 if (w_ls->n_schur_compls > 0) {
 // !!! koncepce Neighbouringu je tak prohnila, ze neni vubec jasne, jestli
 // pro jeden element sousedi max s jednou edge a naopak takze musim spolehat jen na to co je
 // v e_col
 FOR_ELEMENTS( ele ) {
 for(i=0;i<ele->n_sides;i++) w_ls->nonzeros[ele->side[i]->edge->c_row]+=ele->e_row_count;
 for(i=0;i<ele->e_row_count;i++) w_ls->nonzeros[ele->e_col[i]]+=ele->n_sides+ele->e_row_count-1;
 }
 }
 // -C'*A-*B block conect edge with its elements = n_sides nz
 // -C'*A-*C block conect all edges of every element = n_sides*(dim of sides +1)nz (not counting diagonal)
 FOR_EDGES ( edg ) {
 w_ls->nonzeros[row++] += edg->n_sides*(edg->side[0]->dim+1+1);
 }
 }
 }
 */

/*******************************************************************************
 * COMPOSE WATER MH MATRIX WITHOUT SCHUR COMPLEMENT
 ******************************************************************************/

void make_schur0(DarcyFlowMH *w) {
    int i_loc, el_row;
    Element *ele;
    Vec aux;

    Timing *asm_time = timing_create("PREALLOCATION", PETSC_COMM_WORLD);
    if (w->schur0 == NULL) { // create Linear System for MH matrix
        xprintf( Msg, "Allocating MH matrix for water model ... \n " );

        if (w->solver->type == PETSC_MATIS_SOLVER)
            w->schur0 = new LinSys_MATIS(w->lsize, w->ndof_loc,
                    w->global_row_4_sub_row);
        else
            w->schur0 = new LinSys_MPIAIJ(w->lsize);
        w->schur0->set_symmetric();
        w->schur0->start_allocation();
        mh_abstract_assembly(w); // preallocation

    }
    timing_reuse(asm_time, "ASSEMBLY");
    xprintf( Msg, "Assembling MH matrix for water model ... \n " );

    w->schur0->start_add_assembly(); // finish allocation and create matrix
    mh_abstract_assembly(w); // fill matrix
    w->schur0->finalize();
    w->schur0->view_local_matrix();
    timing_destroy(asm_time);

    // add time term

    /*
     for(i_loc=0;i_loc<w->el_ds->lsize;i_loc++ ) {
     ele=w->mesh->element_hash[w->el_4_loc[i_loc]];
     el_row=w->el_row_4_id[ele->id];
     xprintf(Msg,"tdiag: %d %f %f\n",el_row,ele->tAddDiag,ele->tAddRHS);
     LSMatSetValue(w->schur0,el_row,el_row,ele->tAddDiag);
     LSVecSetValue(w->schur0,el_row,ele->tAddRHS);
     }
     */
    //MatView(mtx->mtx,	PETSC_VIEWER_STDOUT_SELF );

    /*
     VecCreateMPI(PETSC_COMM_WORLD,w->lsize,PETSC_DETERMINE,&(aux));
     MatGetDiagonal(w->schur0->A,aux);
     MyVecView(aux,w->old_4_new,"A.dat");
     */

}

//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
void destroy_water_linsys(DarcyFlowMH *w) {
    if (w->schur2 != NULL)
        delete w->schur2;
    if (w->schur1 != NULL)
        delete w->schur1;
    delete w->schur0;

    if (w->solver->type == PETSC_MATIS_SOLVER) {
        xfree(w->global_row_4_sub_row);
    }
    xfree( w );
}

/*******************************************************************************
 * COMPUTE THE FIRST (A-block) SCHUR COMPLEMENT
 ******************************************************************************/
// paralellni verze musi jeste sestrojit index set bloku A, to jde pres:
// lokalni elementy -> lokalni sides -> jejich id -> jejich radky
// TODO: reuse IA a Schurova doplnku
void make_schur1(DarcyFlowMH *w) {
    Mesh* mesh = w->mesh;
    Mat IA;
    ElementFullIter ele = ELEMENT_FULL_ITER(NULL);
    int i_loc, nsides, i, side_rows[4], ierr, el_row;
    double det;
    F_ENTRY;
    Timing *schur1_time = timing_create("Schur 1", PETSC_COMM_WORLD);

    // create Inverse of the A block
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, w->side_ds->lsize(),
            w->side_ds->lsize(), PETSC_DETERMINE, PETSC_DETERMINE, 4,
            PETSC_NULL, 0, PETSC_NULL, &(IA));
    MatSetOption(IA, MAT_SYMMETRIC, PETSC_TRUE);

    for (i_loc = 0; i_loc < w->el_ds->lsize(); i_loc++) {
        ele = w->mesh->element(w->el_4_loc[i_loc]);
        el_row = w->row_4_el[w->el_4_loc[i_loc]];
        nsides = ele->n_sides;
        if (ele->loc_inv == NULL) {
            ele->loc_inv = (double *) malloc(nsides * nsides * sizeof(double));
            det = MatrixInverse(ele->loc, ele->loc_inv, nsides);
            if (fabs(det) < NUM_ZERO) {
                xprintf(Warn,"Singular local matrix of the element %d\n",ele.id());
                PrintSmallMatrix(ele->loc, nsides);
                xprintf(Err,"det: %30.18e \n",det);
            }
        }
        for (i = 0; i < nsides; i++)
            side_rows[i] = w->side_row_4_id[ele->side[i]->id] // side row in MH matrix
                    - w->rows_ds->begin() // local side number
                    + w->side_ds->begin(); // side row in IA matrix
        MatSetValues(IA, nsides, side_rows, nsides, side_rows, ele->loc_inv,
                INSERT_VALUES);
    }

    MatAssemblyBegin(IA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(IA, MAT_FINAL_ASSEMBLY);

    w->schur1 = new SchurComplement(w->schur0, IA);
    w->schur1->form_schur();
    w->schur1->set_spd();

    timing_destroy(schur1_time);
}

/*******************************************************************************
 * COMPUTE THE SECOND (B-block) SCHUR COMPLEMENT
 ******************************************************************************/
void make_schur2(DarcyFlowMH *w) {
    Mat IA;
    Vec Diag, DiagB;
    PetscScalar *vDiag;
    int ierr, loc_el_size;
    F_ENTRY;
    Timing *schur2_time = timing_create("Schur 2", PETSC_COMM_WORLD);
    // create Inverse of the B block ( of the first complement )


    // get subdiagonal of local size == loc num of elements
    loc_el_size = w->el_ds->lsize();
    VecCreateMPI(PETSC_COMM_WORLD, w->schur1->get_system()->vec_lsize(),
            PETSC_DETERMINE, &Diag);
    MatGetDiagonal(w->schur1->get_system()->get_matrix(), Diag); // get whole diagonal
    VecGetArray(Diag,&vDiag);
    // define sub vector of B-block diagonal
    VecCreateMPIWithArray(PETSC_COMM_WORLD, loc_el_size, PETSC_DETERMINE,
            vDiag, &DiagB);
    // compute inverse
    VecReciprocal(DiagB);
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, loc_el_size, loc_el_size,
            PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL, 0, PETSC_NULL,
            &(IA)); // construct matrix
    MatDiagonalSet(IA, DiagB, INSERT_VALUES);
    VecDestroy(DiagB); // clean up
    VecRestoreArray(Diag,&vDiag);
    VecDestroy(Diag);

    w->schur2 = new SchurComplement(w->schur1->get_system(), IA);
    w->schur2->form_schur();
    w->schur2->scale(-1.0);
    w->schur2->set_spd();

    timing_destroy(schur2_time);
}

//-----------------------------------------------------------------------------
// vim: set cindent:
