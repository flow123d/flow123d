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
 * @ingroup flow
 * @brief  Setup and solve linear system of mixed-hybrid discretization of the linear
 * porous media flow with possible preferential flow in fractures and chanels.
 *
 */

#include "petscmat.h"
#include "petscviewer.h"
#include "petscao.h"
#include "petscerror.h"


#include "system/system.hh"
#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "la/distribution.hh"
#include "flow/darcy_flow_mh.hh"
#include "la/linsys.hh"
#include "la/solve.h"
#include "la/schur.hh"
#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"
#include "field_p0.hh"
#include "flow/local_matrix.h"

#include "flow/darcy_flow_mh_output.hh"

#include <limits>
#include <set>
#include <vector>
#include <iostream>
#include <iterator>


Input::Type::AbstractRecord & DarcyFlowMH::get_input_type()
{
    using namespace Input::Type;
    static AbstractRecord rec("DarcyFlowMH", "Mixed-Hybrid  solver for saturated Darcy flow.");

    if (!rec.is_finished()) {
        rec.finish();

        DarcyFlowMH_Steady::get_input_type();

        rec.no_more_descendants();
    }
    return rec;
}





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
DarcyFlowMH_Steady::DarcyFlowMH_Steady(TimeMarks &marks, Mesh &mesh_in, MaterialDatabase &mat_base_in, Input::Record in_rec)
: DarcyFlowMH(marks, mesh_in, mat_base_in)

{
    using namespace Input;
    F_ENTRY;

    int ierr;

    size = mesh_->n_elements() + mesh_->n_sides + mesh_->n_edges();
    n_schur_compls = in_rec.val<int>("n_schurs");
    if ((unsigned int) n_schur_compls > 2) {
        xprintf(Warn,"Invalid number of Schur Complements. Using 2.");
        n_schur_compls = 2;
    }

    solver = new (Solver);
    solver_init(solver, in_rec.val<AbstractRecord>("solver"));

    solution = NULL;
    schur0   = NULL;
    schur1   = NULL;
    schur2   = NULL;
    IA1      = NULL;
    IA2      = NULL;

    Iterator<FilePath> it = in_rec.find<FilePath>("sources_file");
    if (it) {
        sources= new FieldP0<double>(mesh_);
        sources->read_field(*it,string("$Sources"));
    }

    Iterator<string> it_f = in_rec.find<string>("sources_formula");
    if (it) {
        sources= new FieldP0<double>(mesh_);
        sources->setup_from_function(*it_f);
    }

    // time governor
    time_=new TimeGovernor(-1, TimeGovernor::inf_time, *time_marks);

    // init paralel structures
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &(myp));
    ierr += MPI_Comm_size(PETSC_COMM_WORLD, &(np));
    if (ierr)
        xprintf(Err, "Some error in MPI.\n");

    // calculation_mh  - precalculation of some values stored still in mesh_
    {
    struct Side *sde;

    read_boundary(mesh_, in_rec.val<FilePath>("boundary_file") );

    FOR_SIDES(mesh_,sde) {
        calc_side_metrics(sde);
    }

    edge_calculation_mh(mesh_);
    element_calculation_mh(mesh_);
    side_calculation_mh(mesh_);
    boundary_calculation_mh(mesh_);
    local_matrices_mh(mesh_);
    }

    prepare_parallel();

    //side_ds->view();
    //el_ds->view();
    //edge_ds->view();

    make_schur0();

    // prepare Scatter form parallel to sequantial in original numbering
    {
            IS is_par, is_loc;
            int i, si, *loc_idx;
            Edge *edg;

            // create local solution vector
            solution = (double *) xmalloc(size * sizeof(double));
            VecCreateSeqWithArray(PETSC_COMM_SELF, size, solution,
                    &(sol_vec));

            // create seq. IS to scatter par solutin to seq. vec. in original order
            // use essentialy row_4_id arrays
            loc_idx = (int *) xmalloc(size * sizeof(int));
            i = 0;
            FOR_ELEMENTS(mesh_, ele) {
                FOR_ELEMENT_SIDES(ele,si) {
                    loc_idx[i++] = side_row_4_id[ele->side[si]->id];
                }
            }
            FOR_ELEMENTS(mesh_, ele) {
                loc_idx[i++] = row_4_el[ele.index()];
            }
            FOR_EDGES(mesh_, edg) {
                loc_idx[i++] = row_4_edge[edg.index()];
            }
            ASSERT( i==size,"Size of array does not match number of fills.\n");
            //DBGPRINT_INT("loc_idx",size,loc_idx);
            ISCreateGeneral(PETSC_COMM_SELF, size, loc_idx, PETSC_COPY_VALUES, &(is_loc));
            xfree(loc_idx);
            VecScatterCreate(schur0->get_solution(), is_loc, sol_vec,
                    PETSC_NULL, &par_to_all);
            ISDestroy(&(is_loc));
        }
    solution_changed_for_scatter=true;

}

Input::Type::Record & DarcyFlowMH_Steady::get_input_type()
{
    using namespace Input::Type;
    static Record rec("Steady_MH", "Mixed-Hybrid  solver for STEADY saturated Darcy flow.");

    if (!rec.is_finished()) {
        rec.derive_from(DarcyFlowMH::get_input_type());
        rec.declare_key("n_schurs", Integer(0,2), Default("2"),
                "Number of Schur complements to perform when solving MH sytem.");
        rec.declare_key("sources_file", FileName::input(),
                "File with water source field.");
        rec.declare_key("sources_formula", String(),
                "Formula to determine the source field.");
        rec.declare_key("boundary_file", FileName::input(),Default::obligatory(),
                "File with boundary conditions for MH solver.");
        rec.declare_key("solver", Solver::get_input_type(), Default::obligatory(),
                "Linear solver for MH problem.");
        rec.declare_key("output", DarcyFlowMHOutput::get_input_type(), Default::obligatory(),
                "Parameters of output form MH module.");
        rec.finish();
    }
    return rec;
}




//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void DarcyFlowMH_Steady::update_solution() {


    START_TIMER("SOLVING MH SYSTEM");
    F_ENTRY;

    if (time_->is_end()) return;


    time_->next_time();

    xprintf(Msg, "t: %f (Darcy) dt: %f\n",time_->t(), time_->dt());
    modify_system(); // hack for unsteady model

    switch (n_schur_compls) {
    case 0: /* none */
        solve_system(solver, schur0);
        break;
    case 1: /* first schur complement of A block */
        make_schur1();
        //solve_system(solver, schur1->get_system());
        //schur1->resolve();
        schur1->solve(solver);
        break;
    case 2: /* second schur complement of the max. dimension elements in B block */
        make_schur1();
        make_schur2();

        //mat_count_off_proc_values(schur2->get_system()->get_matrix(),schur2->get_system()->get_solution());
        solve_system(solver, schur2->get_system());

        schur2->resolve();
        schur1->resolve();
        break;
    }
    postprocess();

    //int rank;
    //MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    //if ( rank == 0 ) {
    //    PetscViewer solViewer;
    //    PetscViewerASCIIOpen( PETSC_COMM_SELF, "sol.m", &solViewer );
    //    PetscViewerSetFormat(solViewer,PETSC_VIEWER_ASCII_MATLAB);
    //    VecView( sol_vec, solViewer );
    //    PetscViewerDestroy(solViewer);
    //}

    solution_changed_for_scatter=true;
}

void DarcyFlowMH_Steady::postprocess() {
    if (sources == NULL)
        return;

    int i_loc, side_rows[4];
    double values[4];
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    ;

    // modify side fluxes in parallel
    // for every local edge take time term on digonal and add it to the corresponding flux
    for (i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        ele = mesh_->element(el_4_loc[i_loc]);
        FOR_ELEMENT_SIDES(ele,i) {
            side_rows[i] = side_row_4_id[ele->side[i]->id];
            values[i] = -1.0 * ele->volume * sources->element_value(ele.index()) / ele->n_sides;
        }
        VecSetValues(schur0->get_solution(), ele->n_sides, side_rows, values, ADD_VALUES);
    }
    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());
}


void  DarcyFlowMH_Steady::get_solution_vector(double * &vec, unsigned int &vec_size)
{
    // TODO: make class for vectors (wrapper for PETSC or other) derived from LazyDependency
    // and use its mechanism to manage dependency between vectors
    if (solution_changed_for_scatter) {
        // scatter solution to all procs
        VecScatterBegin(par_to_all, schur0->get_solution(), sol_vec,
                INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(par_to_all, schur0->get_solution(), sol_vec,
                INSERT_VALUES, SCATTER_FORWARD);
        solution_changed_for_scatter=false;
    }

    vec=solution;
    vec_size = size;
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}

void  DarcyFlowMH_Steady::get_parallel_solution_vector(Vec &vec)
{
    vec=schur0->get_solution();
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
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
void DarcyFlowMH_Steady::assembly_steady_mh_matrix() {
    LinSys *ls = schur0;
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    int el_row, side_row, edge_row;
    int tmp_rows[100];
    int i, i_loc, nsides, li, si;
    int side_rows[4], edge_rows[4]; // rows for sides and edges of one element
    double f_val;
    double zeros[1000]; // to make space for second schur complement, max. 10 neighbour edges of one el.
    double minus_ones[4] = { -1.0, -1.0, -1.0, -1.0 };
    F_ENTRY;

    //DBGPRINT_INT("side_row_4_id",mesh->max_side_id+1,side_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh->max_elm_id+1,el_row_4_id);
    //DBGPRINT_INT("edge_row_4_id",mesh->max_edg_id+1,edge_row_4_id);
    //DBGPRINT_INT("el_id_4_loc",el_ds->lsize(),el_id_4_loc);

    SET_ARRAY_ZERO(zeros,1000);
    for (i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        ele = mesh_->element(el_4_loc[i_loc]);
        el_row = row_4_el[el_4_loc[i_loc]];
        nsides = ele->n_sides;
        for (i = 0; i < nsides; i++) {
            side_row = side_rows[i] = side_row_4_id[ele->side[i]->id];
            edge_row = edge_rows[i] = row_4_edge[mesh_->edge.index(ele->side[i]->edge)];
            // set block C and C': side-edge, edge-side
            ls->mat_set_value(side_row, edge_row, ele->side[i]->c_val);
            ls->mat_set_value(edge_row, side_row, ele->side[i]->c_val);
        }
        // set block A: side-side on one element - block diagonal matrix
        ls->mat_set_values(nsides, side_rows, nsides, side_rows, ele->loc);
        // set block B, B': element-side, side-element
        ls->mat_set_values(1, &el_row, nsides, side_rows, minus_ones);
        ls->mat_set_values(nsides, side_rows, 1, &el_row, minus_ones);
        // set RHS for sides - dirichlet BC; gravity term ?
        // RHS for elements - neuman BC
        ls->rhs_set_values(nsides, side_rows, ele->rhs);

        // set sources
        if (sources != NULL) {
            ls->rhs_set_value(el_row, -1.0 * ele->volume * sources->element_value(ele.index()));
        }

        // D block: non-compatible conections and diagonal: element-element
        for (i = 0; i < ele->d_row_count; i++)
            tmp_rows[i] = row_4_el[ele->d_el[i]];
        ls->mat_set_values(1, &el_row, ele->d_row_count, tmp_rows, ele->d_val);
        // E',E block: compatible connections: element-edge
        for (i = 0; i < ele->e_row_count; i++)
            tmp_rows[i] = row_4_edge[ele->e_edge_idx[i]];
        ls->mat_set_values(1, &el_row, ele->e_row_count, tmp_rows, ele->e_val);
        ls->mat_set_values(ele->e_row_count, tmp_rows, 1, &el_row, ele->e_val);

        // add virtual values for schur complement allocation
        switch (n_schur_compls) {
        case 2:
            if (ele->d_row_count > 1) {
                xprintf(Warn,"Can not use second Schur complement for problem with non-compatible connections.\n");
                n_schur_compls = 1;
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
    // also Dirichlet BC
    for (i_loc = 0; i_loc < edge_ds->lsize(); i_loc++) {
        edge_row = row_4_edge[edge_4_loc[i_loc]];
        EdgeFullIter edg = mesh_->edge(edge_4_loc[i_loc]);

        //xprintf(Msg,"F: %d %f\n",old_4_new[edge_row],edg->f_val);
        ls->mat_set_value(edge_row, edge_row, edg->f_val);
        ls->rhs_set_value(edge_row, edg->f_rhs);
    }
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

void DarcyFlowMH_Steady::make_schur0() {
    int i_loc, el_row;
    Element *ele;
    Vec aux;

    START_TIMER("PREALLOCATION");

    if (schur0 == NULL) { // create Linear System for MH matrix

        if (solver->type == PETSC_MATIS_SOLVER) 
            schur0 = new LinSys_MATIS( global_row_4_sub_row );
	
        else
            schur0 = new LinSys_MPIAIJ(lsize);
        schur0->set_symmetric();
        schur0->start_allocation();
        assembly_steady_mh_matrix(); // preallocation

    }

    END_TIMER("PREALLOCATION");

    START_TIMER("ASSEMBLY");

    schur0->start_add_assembly(); // finish allocation and create matrix
    assembly_steady_mh_matrix(); // fill matrix
    schur0->finalize();

    //schur0->view_local_matrix();
    //PetscViewer myViewer;
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matis.m",&myViewer);
    //PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
    //MatView( schur0->get_matrix( ), myViewer );
    //PetscViewerDestroy(myViewer);



    // add time term

}

//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
DarcyFlowMH_Steady::~DarcyFlowMH_Steady() {
    if (schur2 != NULL) delete schur2;
    if (schur1 != NULL) delete schur1;

    if ( IA1 != NULL ) MatDestroy( &(IA1) );
    if ( IA2 != NULL ) MatDestroy( &(IA2) );

    delete schur0;
    delete solver;
}

/*******************************************************************************
 * COMPUTE THE FIRST (A-block) SCHUR COMPLEMENT
 ******************************************************************************/
// paralellni verze musi jeste sestrojit index set bloku A, to jde pres:
// lokalni elementy -> lokalni sides -> jejich id -> jejich radky
// TODO: reuse IA a Schurova doplnku
void DarcyFlowMH_Steady::make_schur1() {
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    int i_loc, nsides, i, side_rows[4], ierr, el_row;
    double det;
    PetscErrorCode err;

    F_ENTRY;
    START_TIMER("Schur 1");


    // check type of LinSys
    if      (schur0->type == LinSys::MAT_IS)
    {
       // create mapping for PETSc

       err = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,
               side_ds->lsize(),
               side_id_4_loc, PETSC_COPY_VALUES, &map_side_local_to_global);
       ASSERT(err == 0,"Error in ISLocalToGlobalMappingCreate.");

       err = MatCreateIS(PETSC_COMM_WORLD,  side_ds->lsize(), side_ds->lsize(), side_ds->size(), side_ds->size(), map_side_local_to_global, &IA1);
       ASSERT(err == 0,"Error in MatCreateIS.");

       MatSetOption(IA1, MAT_SYMMETRIC, PETSC_TRUE);

       for (i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
           ele = mesh_->element(el_4_loc[i_loc]);
           el_row = row_4_el[el_4_loc[i_loc]];
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
	   /* print the matrix */
	   //int j;
	   //xprintf(Msg,"Local element inverse: \n ");
           //for (i = 0; i < nsides; i++) {
           //   for (j = 0; j < nsides; j++) 
	   //      xprintf(Msg, " %14.6f ", ele->loc_inv[i*nsides + j]);
	   //   xprintf(Msg, " \n ");
	   //}

           for (i = 0; i < nsides; i++)
               side_rows[i] = ele->side[i]->id; // side ID
                           // - rows_ds->begin(); // local side number
                           // + side_ds->begin(); // side row in IA1 matrix
           MatSetValues(IA1, nsides, side_rows, nsides, side_rows, ele->loc_inv,
                        INSERT_VALUES);
       }
    }
    else if (schur0->type == LinSys::MAT_MPIAIJ)
    {
       if (schur1 == NULL) {
       // create Inverse of the A block
       ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, side_ds->lsize(),
               side_ds->lsize(), PETSC_DETERMINE, PETSC_DETERMINE, 4,
               PETSC_NULL, 0, PETSC_NULL, &(IA1));

       MatSetOption(IA1, MAT_SYMMETRIC, PETSC_TRUE);
       schur1 = new SchurComplement(schur0, IA1);
       }

       for (i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
           ele = mesh_->element(el_4_loc[i_loc]);
           el_row = row_4_el[el_4_loc[i_loc]];
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
               side_rows[i] = side_row_4_id[ele->side[i]->id] // side row in MH matrix
                       - rows_ds->begin() // local side number
                       + side_ds->begin(); // side row in IA1 matrix
           MatSetValues(IA1, nsides, side_rows, nsides, side_rows, ele->loc_inv,
                   INSERT_VALUES);
       }
    }

    MatAssemblyBegin(IA1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(IA1, MAT_FINAL_ASSEMBLY);

    schur1->form_schur();
    schur1->set_spd();
}

/*******************************************************************************
 * COMPUTE THE SECOND (B-block) SCHUR COMPLEMENT
 ******************************************************************************/
void DarcyFlowMH_Steady::make_schur2() {
    PetscScalar *vDiag;
    int ierr, loc_el_size;
    F_ENTRY;
    START_TIMER("Schur 2");
    // create Inverse of the B block ( of the first complement )

    loc_el_size = el_ds->lsize();

    if (schur2 == NULL) {
      // get subdiagonal of local size == loc num of elements
      VecCreateMPI(PETSC_COMM_WORLD, schur1->get_system()->vec_lsize(),
              PETSC_DETERMINE, &diag_schur1);
      ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, loc_el_size, loc_el_size,
              PETSC_DETERMINE, PETSC_DETERMINE, 1, PETSC_NULL, 0, PETSC_NULL,
              &(IA2)); // construct matrix
      VecGetArray(diag_schur1,&vDiag);
      // define sub vector of B-block diagonal
      VecCreateMPIWithArray(PETSC_COMM_WORLD, loc_el_size, PETSC_DETERMINE,
              vDiag, &diag_schur1_b);
      VecRestoreArray(diag_schur1,&vDiag);
      schur2 = new SchurComplement(schur1->get_system(), IA2);

    }

    MatGetDiagonal(schur1->get_system()->get_matrix(), diag_schur1); // get whole diagonal
    // compute inverse of B-block
    VecReciprocal(diag_schur1_b);
    MatDiagonalSet(IA2, diag_schur1_b, INSERT_VALUES);

    schur2->form_schur();
    schur2->scale(-1.0);
    schur2->set_spd();
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

    FOR_EDGES(mesh, edg) {

        // skip non-local edges
        if (!edistr.is_local(edg.index()))
            continue;

        e_weight = edge_dim_weights[edg->side[0]->element->dim - 1];
        // for all connected elements
        FOR_EDGE_SIDES( edg, li ) {
            ASSERT(NONULL(edg->side[li]),"NULL side of edge.");
            ele = edg->side[li]->element;
            ASSERT(NONULL(ele),"NULL element of side.");

            // for sides of connected element, excluding edge itself
            FOR_ELEMENT_SIDES( ele, si ) {
                eid = mesh->edge.index(ele->side[si]->edge);
                if (eid != edg.index())
                    graph->set_edge(edg.index(), eid, e_weight);
            }

            // include connections from lower dim. edge
            // to the higher dimension
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb; i_neigh++) {
                eid = mesh->edge.index(ele->neigh_vb[i_neigh]->edge);
                graph->set_edge(edg.index(), eid, e_weight);
                graph->set_edge(eid, edg.index(), e_weight);
            }
        }
        i_edg++;
    }

    graph->finalize();
}

/**
 * Make connectivity graph of elements of mesh - dual graph: elements vertices of graph.
 */
void make_element_connection_graph(Mesh *mesh, SparseGraph * &graph,bool neigh_on) {

    Distribution edistr=graph->get_distr();

    Edge *edg;
    int li, si, e_idx, i_neigh;
    int i_s, n_s;
    F_ENTRY;

    FOR_ELEMENTS(mesh, ele) {
        //xprintf(Msg,"Element id %d , its index %d.\n",ele.id(), i_ele);

        // skip non-local elements
        if (!edistr.is_local(ele.index()))
            continue;

        // for all connected elements
        FOR_ELEMENT_SIDES( ele, si ) {
            edg = ele->side[si]->edge;

            FOR_EDGE_SIDES( edg, li ) {
                ASSERT(NONULL(edg->side[li]),"NULL side of edge.");
                e_idx = ELEMENT_FULL_ITER(mesh, edg->side[li]->element).index();

                // for elements of connected elements, excluding element itself
                if (e_idx != ele.index()) {
                    graph->set_edge(ele.index(), e_idx);
                }
            }
        }

        // include connections from lower dim. edge
        // to the higher dimension
        if ( neigh_on ) {
            for(i_neigh=0; i_neigh < ele->n_neighs_vb; i_neigh++) {
               n_s = ele->neigh_vb[i_neigh]->edge->n_sides;
               for(i_s=0; i_s < n_s; i_s++) {
                   e_idx=ELEMENT_FULL_ITER(mesh, ele->neigh_vb[i_neigh]->edge->side[i_s]->element).index();
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

    ISCreateGeneral(PETSC_COMM_WORLD, old_ds.lsize(), loc_part, PETSC_COPY_VALUES, &part); // global IS part.
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
    //DBGMSG("Creating global new->old mapping ...\n");
    AOCreateBasicIS(new_numbering, PETSC_NULL, &new_old_ao); // app ordering= new; petsc ordering = old
    for (i = 0; i < size; i++)
        old_4_new[i] = i;
    AOApplicationToPetsc(new_old_ao, size, old_4_new);
    AODestroy(&(new_old_ao));

    // compute id_4_loc
    //DBGMSG("Creating loc.number -> id mapping ...\n");
    i_loc = 0;
    //DBGPRINT_INT("id_4_old",old_ds.lsize(),id_4_old);
    //DBGPRINT_INT("old_4_new",new_ds->lsize(),old_4_new)

    for (i_new = new_ds->begin(); i_new < new_ds->end(); i_new++) {
        //printf("i_new: %d old: %d id: %d i_loc: %d \n",i_new,old_4_new[i_new],i_loc);
        id_4_loc[i_loc++] = id_4_old[old_4_new[i_new]];
    }
    // compute row_4_id
    //DBGMSG("Creating id -> stiffness mtx. row mapping ...\n");
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
void DarcyFlowMH_Steady::make_row_numberings() {
    int i, shift;
    int np = edge_ds->np();
    int edge_shift[np], el_shift[np], side_shift[np];
    unsigned int rows_starts[np];

    int edge_n_id = mesh_->n_edges(),
            el_n_id = mesh_->element.size(),
            side_n_id = mesh_->n_sides;

    // compute shifts on every proc
    shift = 0; // in this var. we count new starts of arrays chunks
    for (i = 0; i < np; i++) {
        side_shift[i] = shift - (side_ds->begin(i)); // subtract actual start of the chunk
        shift += side_ds->lsize(i);
        el_shift[i] = shift - (el_ds->begin(i));
        shift += el_ds->lsize(i);
        edge_shift[i] = shift - (edge_ds->begin(i));
        shift += edge_ds->lsize(i);
        rows_starts[i] = shift;
    }
    //DBGPRINT_INT("side_shift",np,side_shift);
    //DBGPRINT_INT("el_shift",np,el_shift);
    //DBGPRINT_INT("edge_shift",np,edge_shift);
    // apply shifts
    for (i = 0; i < side_n_id; i++) {
        int &what = side_row_4_id[i];
        if (what >= 0)
            what += side_shift[side_ds->get_proc(what)];
    }
    for (i = 0; i < el_n_id; i++) {
        int &what = row_4_el[i];
        if (what >= 0)
            what += el_shift[el_ds->get_proc(what)];

    }
    for (i = 0; i < edge_n_id; i++) {
        int &what = row_4_edge[i];
        if (what >= 0)
            what += edge_shift[edge_ds->get_proc(what)];
    }
    // make distribution of rows
    for (i = np - 1; i > 0; i--)
        rows_starts[i] -= rows_starts[i - 1];


    rows_ds = boost::make_shared<Distribution>(&(rows_starts[0]));
}

// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void DarcyFlowMH_Steady::prepare_parallel() {

    int *loc_part; // optimal (edge,el) partitioning (local chunk)
    int *id_4_old; // map from old idx to ids (edge,el)
    // auxiliary
    Edge *edg;
    Element *el;
    Side *side;
    int i, loc_i;

    int i_neigh;
    int e_idx;
    int i_loc, el_row, side_row, edge_row, nsides;

    PetscErrorCode err;
    F_ENTRY;
    MPI_Barrier(PETSC_COMM_WORLD);

    if (solver->type == PETSC_MATIS_SOLVER) {
        xprintf(Msg,"Compute optimal partitioning of elements.\n");

        // prepare dual graph
        Distribution init_el_ds(Distribution::Block, mesh_->n_elements());  // initial distr.
        SparseGraph *element_graph= new SparseGraphPETSC(init_el_ds);
        loc_part = new int[init_el_ds.lsize()];                                     // partitionig in initial distribution

        make_element_connection_graph(mesh_, element_graph);

	//element_graph->view();

        WARN_ASSERT(element_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

        element_graph->partition(loc_part);
        //DBGPRINT_INT("loc_part",init_el_ds.lsize(),loc_part);

        // prepare parallel distribution of dofs linked to elements
        id_4_old = new int[mesh_->n_elements()];
        i = 0;
        FOR_ELEMENTS(mesh_, el)
            id_4_old[i++] = el.index();
        id_maps(mesh_->element.size(), id_4_old, init_el_ds, loc_part, el_ds,
                el_4_loc, row_4_el);
        //DBGPRINT_INT("el_4_loc",el_ds->lsize(),el_4_loc);
        //xprintf(Msg,"Number of elements in subdomain %d \n",el_ds->lsize());
        delete[] loc_part;
        delete[] id_4_old;

	//el_ds->view();
	//
        //DBGMSG("Compute appropriate edge partitioning ...\n");
        //optimal element part; loc. els. id-> new el. numbering
        Distribution init_edge_ds(Distribution::Localized, mesh_->n_edges());
        // partitioning of edges, edge belongs to the proc of his first element
        // this is not optimal but simple
        loc_part = new int[init_edge_ds.lsize()];
        id_4_old = new int[mesh_->n_edges()];
        {
            loc_i = 0;
            FOR_EDGES(mesh_, edg ) {
                // partition
                e_idx = mesh_->element.index(edg->side[0]->element);
                //xprintf(Msg,"Index of edge: %d first element: %d \n",edgid,e_idx);
                if (init_edge_ds.is_local(edg.index())) {
                    // find (new) proc of the first element of the edge
                    loc_part[loc_i++] = el_ds->get_proc(row_4_el[e_idx]);
                }
                // id array
                id_4_old[edg.index()] = edg.index();
            }
        }
        //    // make trivial part
        //    for(loc_i=0;loc_i<init_el_ds->lsize;loc_i++) loc_part[loc_i]=init_el_ds->myp;
        //DBGPRINT_INT("loc_part",init_edge_ds.lsize(),loc_part);

        id_maps(mesh_->n_edges(), id_4_old, init_edge_ds, loc_part, edge_ds,
                edge_4_loc, row_4_edge);
        delete[] loc_part;
        delete[] id_4_old;


    } else {
        xprintf(Msg,"Compute optimal partitioning of edges.\n");

        SparseGraph *edge_graph = new SparseGraphMETIS(mesh_->n_edges());                     // graph for partitioning
        Distribution init_edge_ds = edge_graph->get_distr();  // initial distr.
        loc_part = new int[init_edge_ds.lsize()];                                     // partitionig in initial distribution

        make_edge_conection_graph(mesh_, edge_graph);
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
        id_4_old = new int[mesh_->n_edges()];
        i = 0;
        FOR_EDGES(mesh_, edg)
            id_4_old[i++] = edg.index();
        id_maps(mesh_->n_edges(), id_4_old, init_edge_ds, (int *) loc_part,
                edge_ds, edge_4_loc, row_4_edge);

        delete[] loc_part;
        delete[] id_4_old;

        //DBGMSG("Compute appropriate element partitioning ...\n");
        //optimal element part; loc. els. id-> new el. numbering
        Distribution init_el_ds(Distribution::Block, mesh_->n_elements());
        // partitioning of elements, element belongs to the proc of his first edge
        // this is not optimal but simple
        loc_part = new int[init_el_ds.lsize()];
        id_4_old = new int[mesh_->n_elements()];
        {
            int i_edg;
            loc_i = 0;
            FOR_ELEMENTS(mesh_,  el ) {
                // partition
                if (init_el_ds.is_local( el.index() )) {
                    // find (new) proc of the first edge of element
                    //DBGMSG("%d %d %d %d\n",iel,loc_i,el->side[0]->edge->id,edge_row_4_id[el->side[0]->edge->id]);
                    i_edg=mesh_->edge.index(el->side[0]->edge); // global index in old numbering
                    loc_part[loc_i++] = edge_ds->get_proc(row_4_edge[i_edg]);

                }
                // id array
                id_4_old[el.index()] = el.index();
            }
        }
        //    // make trivial part
        //    for(loc_i=0;loc_i<init_el_ds->lsize;loc_i++) loc_part[loc_i]=init_el_ds->myp;
        id_maps(mesh_->element.size(), id_4_old, init_el_ds, loc_part, el_ds,
                el_4_loc, row_4_el);
        delete[] loc_part;
        delete[] id_4_old;
    }

    //DBGMSG("Compute side partitioning ...\n");
    //optimal side part; loc. sides; id-> new side numbering
    Distribution init_side_ds(Distribution::Block, mesh_->n_sides);
    // partitioning of sides follows elements
    loc_part = new int[init_side_ds.lsize()];
    id_4_old = new int[mesh_->n_sides];
    {
        int is = 0, iel;
        loc_i = 0;
        FOR_SIDES(mesh_, side ) {
            // partition
            if (init_side_ds.is_local(is)) {
                // find (new) proc of the element of the side
                loc_part[loc_i++] = el_ds->get_proc(
                        row_4_el[mesh_->element.index(side->element)]);
            }
            // id array
            id_4_old[is++] = side->id;
        }
    }
    // make trivial part
    //for(loc_i=0;loc_i<init_side_ds->lsize;loc_i++) loc_part[loc_i]=init_side_ds->myp;

    id_maps(mesh_->n_sides, id_4_old, init_side_ds, loc_part, side_ds,
            side_id_4_loc, side_row_4_id);
    delete [] loc_part;
    delete [] id_4_old;

    /*
     DBGPRINT_INT("edge_id_4_loc",edge_ds->lsize,edge_id_4_loc);
     DBGPRINT_INT("el_4_loc",el_ds->lsize,el_4_loc);
     DBGPRINT_INT("side_id_4_loc",side_ds->lsize,side_id_4_loc);
     DBGPRINT_INT("edge_row_4_id",mesh_->n_edges,edge_row_4_id);
     DBGPRINT_INT("el_row_4_id",mesh_->max_elm_id+1,el_row_4_id);
     DBGPRINT_INT("side_row_4_id",mesh_->max_side_id+1,side_row_4_id);
     */
    // convert row_4_id arrays from separate numberings to global numbering of rows
    //MPI_Barrier(PETSC_COMM_WORLD);
    //DBGMSG("Finishing row_4_id\n");
    //MPI_Barrier(PETSC_COMM_WORLD);
    make_row_numberings();
    //DBGPRINT_INT("edge_row_4_id",mesh_->n_edges,edge_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh_->max_elm_id+1,el_row_4_id);
    //DBGPRINT_INT("side_row_4_id",mesh_->max_side_id+1,side_row_4_id);

    lsize = side_ds->lsize() + el_ds->lsize() + edge_ds->lsize();

    /*
    // make old_4_new
    old_4_new = (int *) malloc((mesh_->n_edges() + mesh_->n_sides
            + mesh_->n_elements()) * sizeof(int));
    i = 0;
    FOR_SIDES( side )
        old_4_new[side_row_4_id[side->id]] = i++;
    FOR_ELEMENTS( el )
        old_4_new[row_4_el[el.index()]] = i++;
    FOR_EDGES(edg)
        old_4_new[edge_row_4_id[edg->id]] = i++;
    */
    // prepare global_row_4_sub_row
    if (solver->type == PETSC_MATIS_SOLVER) {
        //xprintf(Msg,"Compute mapping of local subdomain rows to global rows.\n");

        global_row_4_sub_row = boost::make_shared<LocalToGlobalMap>(rows_ds);

        //
        // ordering of dofs
        // for each subdomain:
        // | velocities (at sides) | pressures (at elements) | L. mult. (at edges) |
        for (i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
            el = mesh_->element(el_4_loc[i_loc]);
            el_row = row_4_el[el_4_loc[i_loc]];

            global_row_4_sub_row->insert( el_row );

            nsides = el->n_sides;
            for (i = 0; i < nsides; i++) {
                side_row = side_row_4_id[el->side[i]->id];
                Edge *edg=el->side[i]->edge; 
		        edge_row = row_4_edge[mesh_->edge.index(edg)];

		        global_row_4_sub_row->insert( side_row );
		        global_row_4_sub_row->insert( edge_row );

		// edge neighbouring overlap
		//if (edg->neigh_vb != NULL) {
		//	int neigh_el_row=row_4_el[mesh_->element.index(edg->neigh_vb->element[0])];
                //      localDofSet.insert( neigh_el_row );
		//}	
            }

            for (i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
                // mark this edge
                edge_row = row_4_edge[mesh_->edge.index(el->neigh_vb[i_neigh]->edge)];
                global_row_4_sub_row->insert( edge_row );
            }
        }
        global_row_4_sub_row->finalize();
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
    //printf("[%d] rows: %d off_rows: %d on: %d off: %d\n",distr.myp(),last-first,n_off_rows,n_on,n_off);
}





// ========================
// unsteady

DarcyFlowMH_Unsteady::DarcyFlowMH_Unsteady(TimeMarks &marks,Mesh &mesh_in, MaterialDatabase &mat_base_in, Input::Record in_rec)
    : DarcyFlowMH_Steady(marks,mesh_in, mat_base_in, in_rec)
{
    delete time_; // delete steady TG

    // time governor
    time_=new TimeGovernor(
            0.0,
            OptGetDbl("Global", "Stop_time", "1.0"),
            *time_marks
            );

    time_->set_permanent_constrain(
            OptGetDbl("Global", "Time_step", "1.0"),
            OptGetDbl("Global", "Time_step", "1.0")
            );
    time_->fix_dt_until_mark();
    setup_time_term();


}

void DarcyFlowMH_Unsteady::setup_time_term() {

    // have created full steady linear system
    // save diagonal of steady matrix
    VecCreateMPI(PETSC_COMM_WORLD,rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
    MatGetDiagonal(schur0->get_matrix(), steady_diagonal);

    // read inital condition

    string file_name=IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Initial", "\\" ));
    INPUT_CHECK( file_name != "\\","Undefined filename with initial pressure.\n");
    VecZeroEntries(schur0->get_solution());

    FieldP0<double> *initial_pressure = new FieldP0<double>(mesh_);
    initial_pressure->read_field("input/pressure_initial.in",string("$PressureHead"));
    double *local_sol=schur0->get_solution_array();

    PetscScalar *local_diagonal;
    VecDuplicate(steady_diagonal,& new_diagonal);
    VecZeroEntries(new_diagonal);
    VecGetArray(new_diagonal,& local_diagonal);

    // apply initial condition and modify matrix diagonal
    // cycle over local element rows
    int i_loc_row, i_loc_el;
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    DBGMSG("Setup with dt: %f\n",time_->dt());
    for (i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
        ele = mesh_->element(el_4_loc[i_loc_el]);
        i_loc_row=i_loc_el+side_ds->lsize();

        // set initial condition
        local_sol[i_loc_row]=initial_pressure->element_value(ele.index());
        // set new diagonal
        local_diagonal[i_loc_row]=-ele->material->stor*ele->measure/time_->dt();
    }
    VecRestoreArray(new_diagonal,& local_diagonal);
    MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);
    delete initial_pressure;

    // set previous solution as copy of initial condition
    VecDuplicate(schur0->get_solution(), &previous_solution);
    VecCopy(schur0->get_solution(), previous_solution);

    // save RHS
    VecDuplicate(schur0->get_rhs(), &steady_rhs);
    VecCopy(schur0->get_rhs(),steady_rhs);

    solution_changed_for_scatter=true;
}

void DarcyFlowMH_Unsteady::modify_system()
{

  if (time_->is_changed_dt()) {
      MatDiagonalSet(schur0->get_matrix(),steady_diagonal, INSERT_VALUES);

      VecScale(new_diagonal, time_->last_dt()/time_->dt());
      MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);
  }

  // modify RHS - add previous solution
  VecPointwiseMult(schur0->get_rhs(),new_diagonal,schur0->get_solution());
  VecAXPY(schur0->get_rhs(),1.0,steady_rhs);

  // swap solutions
  VecSwap(previous_solution, schur0->get_solution());
}


// ========================
// unsteady

DarcyFlowLMH_Unsteady::DarcyFlowLMH_Unsteady(TimeMarks &marks,Mesh &mesh_in, MaterialDatabase &mat_base_in, Input::Record in_rec)
    : DarcyFlowMH_Steady(marks,mesh_in, mat_base_in, in_rec)
{
    delete time_; // delete steady TG
    // time governor
    time_=new TimeGovernor(
            0.0,
            OptGetDbl("Global", "Stop_time", "1.0"),
            *time_marks
            );

    time_->set_permanent_constrain(
            OptGetDbl("Global", "Time_step", "1.0"),
            OptGetDbl("Global", "Time_step", "1.0")
            );
    time_->fix_dt_until_mark();
    setup_time_term();

}
void DarcyFlowLMH_Unsteady::setup_time_term()
{
    // have created full steady linear system
     // save diagonal of steady matrix
     VecCreateMPI(PETSC_COMM_WORLD,rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
     MatGetDiagonal(schur0->get_matrix(), steady_diagonal);

     // read inital condition

     string file_name=IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Initial", "\\" ));
     INPUT_CHECK( file_name != "\\","Undefined filename with initial pressure.\n");
     VecZeroEntries(schur0->get_solution());
     FieldP0<double> *initial_pressure = new FieldP0<double>(mesh_);
     initial_pressure->read_field("input/pressure_initial.in",string("$PressureHead"));

     VecDuplicate(steady_diagonal,& new_diagonal);

     // apply initial condition and modify matrix diagonal
     // cycle over local element rows
     int i_loc_row, i_loc_el, edge_row;
     ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
     double init_value;

     for (i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
         ele = mesh_->element(el_4_loc[i_loc_el]);
         i_loc_row=i_loc_el+side_ds->lsize();

         init_value=initial_pressure->element_value(ele.index());

         FOR_ELEMENT_SIDES(ele,i) {
             edge_row = row_4_edge[mesh_->edge.index(ele->side[i]->edge)];
             // set new diagonal
             VecSetValue(new_diagonal,edge_row,-ele->material->stor*ele->volume /time_->dt()/ele->n_sides,ADD_VALUES);
             // set initial condition
             VecSetValue(schur0->get_solution(),edge_row,init_value/ele->n_sides,ADD_VALUES);
         }
     }
     VecAssemblyBegin(new_diagonal);
     VecAssemblyBegin(schur0->get_solution());
     VecAssemblyEnd(new_diagonal);
     VecAssemblyEnd(schur0->get_solution());

     delete initial_pressure;
     MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);

     // set previous solution as copy of initial condition
     VecDuplicate(schur0->get_solution(), &previous_solution);
     VecCopy(schur0->get_solution(), previous_solution);

     // save RHS
     VecDuplicate(schur0->get_rhs(), &steady_rhs);
     VecCopy(schur0->get_rhs(),steady_rhs);

     // auxiliary vector for time term
     VecDuplicate(schur0->get_rhs(), &time_term);

     solution_changed_for_scatter=true;

}

void DarcyFlowLMH_Unsteady::modify_system()
{
    if (time_->is_changed_dt()) {
        MatDiagonalSet(schur0->get_matrix(),steady_diagonal, INSERT_VALUES);

        VecScale(new_diagonal, time_->last_dt()/time_->dt());
        MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);
    }

  // modify RHS - add previous solution
  VecPointwiseMult(schur0->get_rhs(),new_diagonal,schur0->get_solution());
  VecAXPY(schur0->get_rhs(),1.0,steady_rhs);

  // swap solutions
  VecSwap(previous_solution, schur0->get_solution());
}

// TODO: make this operating on parallel solution
// i.e. access from elements to edge values (possibly by constructing specific matrix)

// is it really necessary what is natural value of element pressures ?
// Since
void DarcyFlowLMH_Unsteady::postprocess()
{
  int i_loc,side_row,loc_edge_row,i;
  Edge* edg;
  ElementIter ele;
  double new_pressure, old_pressure,time_coef;

  PetscScalar *loc_prev_sol;
  VecGetArray(previous_solution, &loc_prev_sol);

  // modify side fluxes in parallel
  // for every local edge take time term on digonal and add it to the corresponding flux
  for (i_loc =0 ; i_loc < edge_ds->lsize(); i_loc++) {

      EdgeFullIter edg = mesh_->edge(edge_4_loc[i_loc]);
      loc_edge_row=side_ds->lsize()+el_ds->lsize()+i_loc;

      new_pressure=(schur0->get_solution_array())[loc_edge_row];
      old_pressure=loc_prev_sol[loc_edge_row];
      FOR_EDGE_SIDES(edg,i) {
          ele=edg->side[i]->element;
          side_row=side_row_4_id[edg->side[i]->id];
          time_coef=-ele->material->stor*ele->volume /time_->dt()/ele->n_sides;
          VecSetValue(schur0->get_solution(),side_row,time_coef*(new_pressure-old_pressure),ADD_VALUES);
      }
  }
  VecRestoreArray(previous_solution, &loc_prev_sol);

  VecAssemblyBegin(schur0->get_solution());
  VecAssemblyEnd(schur0->get_solution());

  DarcyFlowMH_Steady::postprocess();
}

//-----------------------------------------------------------------------------
// vim: set cindent:
