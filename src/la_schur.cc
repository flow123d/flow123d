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
 * @file
 * @brief  Assembly explicit Schur complement for the given linear system.
 * Provides method for resolution of the full original vector of unknowns.
 *
 * Aim: Explicit schur should be faster then implicit, i.e.
 *
 * @todo
 * - vyresit navaznost na lin sys - solve a export seq vektoru, redukce ... ?
 * - inv_a - predava se pri konstrukci, ale neumoznuje jeji reuse - aktualizaci assemblace
 *   resp. nutno si na ni drzet ukazatel venku ... reseni ?
 * - ? remove old_4_new - just for LSView
 * - automatic preallocation
 * - eliminated block given by IS
 * - in place Schur
 * - ? nemodifikovat puvodni system, leda skrze jeho metody
 */

#include <petscvec.h>
#include <petscmat.h>

#include <par_distribution.hh>
#include <solve.h>
#include <system.hh>
#include <la_linsys.hh>
#include <la_schur.hh>

/**
 *  Create Schur complement system.
 *  @param[in] Orig  : original system
 *  @param[in] inv_a : inversion of the A block
 *  @param[in] ia    : index set of the A block,
 *   default continuous given by inv_a:
 *   proc    1      2      3
 *
 *   Orig:   ****** ****** ****
 *   IA  :   ***    **     ***
 *
 */


SchurComplement :: SchurComplement(LinSys *Orig, Mat inv_a, IS ia)
: Orig(Orig), IA(inv_a), IsA(ia)
{
    PetscScalar *rhs_array, *sol_array;
    int i, m,n;
    IS ISBloc;

    int orig_first,orig_sub_size, orig_lsize;

    // MATIS vars
    PetscErrorCode err;
    Mat orig_mat_sub;
    int *global_row_4_sub_row_new;
    int SizeA, locSizeB_vec;

    F_ENTRY;

    // check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {
       // it is assumed that Schur complement may be formed block-wise

       // get local submatrix of A_inverse
       err = MatISGetLocalMat(IA, &IA_sub);
       ASSERT(err == 0,"Error in MatISGetLocalMat.");

       // get local submatrix of original system
       err = MatISGetLocalMat(Orig->get_matrix(), &orig_mat_sub);
       ASSERT(err == 0,"Error in MatISGetLocalMat.");

       // find original size
       err = MatGetSize(orig_mat_sub, &m, &n);
       ASSERT(m == n,"Assumed square matrix.");

       // set size
       orig_sub_size = m;
       DBGMSG("orig block size %d",orig_sub_size);

       // find inverted block size
       err = MatGetSize(IA_sub, &m, &n);
       ASSERT(m == n,"Assumed square matrix.");

       // set size
       locSizeA = m;
       DBGMSG("A block size %d",locSizeA);

       orig_first = 0;

       // create A block index set
       if (IsA == NULL) {
           // assume 'a_inv->local_size' be local part of A block
           ISCreateStride(PETSC_COMM_SELF,locSizeA,orig_first,1,&IsA);
       } 
       fullIsA = NULL;
       // create B block index set
       locSizeB = orig_sub_size-locSizeA;

       orig_lsize = Orig->vec_lsize();
       locSizeB_vec = orig_lsize - locSizeA;
       ISCreateStride(PETSC_COMM_WORLD,locSizeB,orig_first+locSizeA,1,&IsB);
       fullIsB = NULL;

    }
    else if (Orig->type == LinSys::MAT_MPIAIJ)
    {
       // get distribution of original marix
       MatGetOwnershipRange(Orig->get_matrix(),&orig_first,PETSC_NULL);
       MatGetLocalSize(Orig->get_matrix(),&orig_lsize,PETSC_NULL);

       // create A block index set
       if (IsA == NULL) {
           // assume 'a_inv->local_size' be local part of A block
           MatGetLocalSize(IA,&locSizeA,PETSC_NULL);
           ISCreateStride(PETSC_COMM_WORLD,locSizeA,orig_first,1,&IsA);
       } else {
           ISGetLocalSize(IsA, &locSizeA);
       }
       ISAllGather(IsA,&fullIsA);

       // create B block index set
       locSizeB = orig_lsize-locSizeA;
       ISCreateStride(PETSC_COMM_WORLD,locSizeB,orig_first+locSizeA,1,&IsB);
       ISAllGather(IsB,&fullIsB);
    }


    // indicate first construction
    mat_reuse=MAT_INITIAL_MATRIX;;

    //DBGMSG("ISs :\n");
    //ISView(Schur->IsA,PETSC_VIEWER_STDOUT_WORLD);
    //ISView(Schur->IsB,PETSC_VIEWER_STDOUT_WORLD);
    //ISView(Schur->fullIsA,PETSC_VIEWER_STDOUT_SELF);
    //ISView(Schur->fullIsB,PETSC_VIEWER_STDOUT_SELF);


    // create complement system
    // TODO: introduce LS as true object, clarify its internal states
    // create RHS sub vecs RHS1, RHS2
    VecGetArray(Orig->get_rhs(),&rhs_array);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeA,PETSC_DETERMINE,rhs_array,&(RHS1));

    // create Solution sub vecs Sol1, Compl->solution
    VecGetArray(Orig->get_solution(),&sol_array);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeA,PETSC_DETERMINE,sol_array,&(Sol1));

    if      (Orig->type == LinSys::MAT_IS)
    {
       VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB_vec,PETSC_DETERMINE,rhs_array+locSizeA,&(RHS2));
       VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB_vec,PETSC_DETERMINE,sol_array+locSizeA,&(Sol2));
    }
    else if (Orig->type == LinSys::MAT_MPIAIJ)
    {
       VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB,PETSC_DETERMINE,rhs_array+locSizeA,&(RHS2));
       VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB,PETSC_DETERMINE,sol_array+locSizeA,&(Sol2));
    }
    VecRestoreArray(Orig->get_rhs(),PETSC_NULL);


    VecRestoreArray(Orig->get_solution(),PETSC_NULL);


    // check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {
       // find size of block A
       err = MatGetSize(inv_a, &m, &n);
       ASSERT(m == n,"Assumed square matrix.");
       SizeA = m;

       global_row_4_sub_row_new = new int[locSizeB];
       for (i = 0; i < locSizeB; i++)
       { 
	  global_row_4_sub_row_new[i] = Orig->subdomain_indices[locSizeA + i] - SizeA;
       }

       Compl = new LinSys_MATIS(locSizeB_vec, locSizeB, global_row_4_sub_row_new, sol_array+locSizeA);

       delete[] global_row_4_sub_row_new;
    }
    else if (Orig->type == LinSys::MAT_MPIAIJ)
    {
       Compl = new LinSys_MPIAIJ(locSizeB, sol_array+locSizeA);
    }



    // TODO: have old_4_new as a mapping inicialized by onother one and
    // parallel shift, use it only if there is not NODEBUG
    // set old_4_new
    /*
    Compl->old_4_new=(int*)malloc(Compl->size*sizeof(int));
    for(i=0;i<Compl->size;i++) {
        n=Compl->ds->get_proc(i);
        //xprintf(Msg,"i s ls:%d %d %d\n",i,Orig->ds->starts[n+1],DS_LSIZE(Schur->Compl->ds,n));
        Compl->old_4_new[i]=Orig->old_4_new[i+Orig->ds->begin(n+1)-Compl->ds->begin(n+1)];
    }
    */
    /*if (Orig->ds->myp==0) {
            DBGMSG("Make Schur compl old_4_new from:\n");
            DBGPRINT_INT("orig old_4_new",Orig->size,Orig->old_4_new);
            DBGPRINT_INT("new old_4_new",Schur->Compl->size,Schur->Compl->old_4_new);
        }
        MPI_Barrier(PETSC_COMM_WORLD);*/

    state=created;

}

/**
 * @brief Form Schur complement. Call solve. Resolve original solution.
 * TODO: should be better when LinSys is full object.
 *
 */

void SchurComplement::solve(Solver *solver)
{
    if (state == created) form_schur();
    ASSERT(state == formed, "Object in wrong state.\n" );

    solve_system(solver, get_system() );
    if      (Orig->type == LinSys::MAT_IS)
    {
        xprintf( MsgDbg, "We deal with LinSys based on MATIS matrix ... not supported for Schur \n " );
        MPI_Barrier(PETSC_COMM_WORLD);
        fflush(stdout);
        MPI_Barrier(PETSC_COMM_WORLD);
        exit(EXIT_FAILURE);
    }
    resolve();

}

void SchurComplement::set_spd()
{Compl->set_positive_definite();}

/**
 *  COMPUTE A SCHUR COMPLEMENT OF A PETSC MATRIX
 *
 *  given symmetric original matrix Orig has form
 *  A  B     x_1    RHS_1
 *  B' C  *  x_2  = RHS_2
 *  where the first block is given by index set IsA, and the second block by IsB
 *  user has to provide inverse IA of the A-block
 *  we suppose that original matrix have non-zero pattern for the schur complement
 *
 *  we return: Shur - schur complement, ShurRHS - RHS of the complemented system:
 *  (B' * IA * B - C) * x_2 = (B' * IA * RHS_1 - RHS_2)
 *  IAB - a matrix to compute eliminated part of the solution:
 *  x_1 = IA * RHS_1 - IAB * x_2
 *
 *  Actually as B' is stored separetly, the routine can be used also for nonsymetric original
 * system
 *
 */

void SchurComplement::form_schur()
{
    // MATIS vars
    PetscScalar *rhs_array_old, *rhs_array_new, *rhs2_array_new;
    PetscErrorCode err;
    Mat orig_mat_sub;
    Vec rhs1_vec;
    Vec rhs_new_vec;
    PetscInt m,n;

    ASSERT(state == created, "Object in wrong state.\n");

    // check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {
       // get local submatrix of original system
       err = MatISGetLocalMat(Orig->get_matrix(), &orig_mat_sub);
       ASSERT(err == 0,"Error in MatISGetLocalMat.");

       // B
       MatGetSubMatrix(orig_mat_sub, IsA, IsB, locSizeB, mat_reuse, &B_sub);

       // A^-1 * B
       MatMatMult(IA_sub, B_sub, mat_reuse, 1.0 ,&(IAB_sub)); // 6/7 - fill estimate

       // B^T
       MatGetSubMatrix(orig_mat_sub, IsB, IsA, locSizeA, mat_reuse, &Bt_sub);

       // B^T*A^-1*B
       MatMatMult(Bt_sub, IAB_sub, mat_reuse, 1.9 ,&(xA_sub)); // 1.1 - fill estimate (PETSC report values over 1.8)

       // get C block
       MatGetSubMatrix(orig_mat_sub, IsB, IsB, locSizeB, mat_reuse, &(Compl->local_matrix));

       // compute complement = (-1)cA+xA = Bt*IA*B - C
       MatScale(Compl->local_matrix,-1.0);
       MatAXPY(Compl->local_matrix, 1, xA_sub, DIFFERENT_NONZERO_PATTERN);
       MatView(Compl->local_matrix,PETSC_VIEWER_STDOUT_SELF);

       // compute the SchurRHS
       //VecGetArray(RHS1,&rhs_array_old);
       //VecCreateSeqWithArray(locSizeA,rhs_array_old,&rhs1_vec)
       //VecRestoreArray(RHS1,PETSC_NULL);


       //VecGetArray(Compl->get_rhs(),&rhs_array_new);
       // find original size
       //err = MatGetSize(IA_sub, &m, &n);
       //ASSERT(m == n,"Assumed square matrix.");
       //VecCreateSeqWithArray(PETSC_COMM_SELF,m,rhs_array_new,&rhs_new_vec);
       //VecRestoreArray(Compl->get_rhs(),PETSC_NULL);
       //MatMultTranspose(IAB_sub,rhs_new_vec,rhs_new_vec);

       //VecGetArray(RHS2,&rhs2_array_new);
       //VecAXPY(rhs_array_new,-1,rhs2_array_new);

    }
    else if      (Orig->type == LinSys::MAT_MPIAIJ)
    {

       DBGMSG("Compute Schur complement of\n");
       //MatView(Schur->Orig->A,PETSC_VIEWER_STDOUT_WORLD);
       DBGMSG("inverse IA:\n");
       //MatView(Schur->IA,PETSC_VIEWER_STDOUT_WORLD);
       // compose Schur complement
       // Petsc need some fill estimate for results of multiplication in form nnz(A*B)/(nnz(A)+nnz(B))
       // for the first Schur compl: IA*B is bounded by ( d*(d+1) )/( d*d+2*d ) <= 5/6 for d<=4
       //                            B'*IA*B bounded by ( (d+1)*(d+1) )/ ( d*(d+1) + d ) ~ 1
       // for the second Schur :      IA*B have fill ratio ~ 1.
       //                            B'*IA*B  ...         ( N/2 *(2*N-1) )/( 2 + 2*N ) <= 1.4
       // nevertheless Petsc does not allows fill ratio below 1. so we use 1.1 for the first
       // and 1.5 for the second multiplication

       // compute IAB=IA*B
       MatGetSubMatrix(Orig->get_matrix(), IsA, fullIsB, locSizeB, mat_reuse, &B);
       DBGMSG(" B:\n");
       //MatView(Schur->B,PETSC_VIEWER_STDOUT_WORLD);
       MatMatMult(IA, B, mat_reuse, 1.0 ,&(IAB)); // 6/7 - fill estimate
       DBGMSG(" IAB:\n");
       //MatView(Schur->IAB,PETSC_VIEWER_STDOUT_WORLD);
       // compute xA=Bt* IAB = Bt * IA * B
       MatGetSubMatrix(Orig->get_matrix(), IsB, fullIsA, locSizeA, mat_reuse, &(Bt));
       MatMatMult(Bt, IAB, mat_reuse, 1.9 ,&(xA)); // 1.1 - fill estimate (PETSC report values over 1.8)
       DBGMSG("xA:\n");
       //MatView(Schur->xA,PETSC_VIEWER_STDOUT_WORLD);

       // get C block
       MatGetSubMatrix(Orig->get_matrix(), IsB, fullIsB, locSizeB, mat_reuse, &(Compl->matrix));
       // compute complement = (-1)cA+xA = Bt*IA*B - C
       MatScale(Compl->get_matrix(),-1.0);
       DBGMSG("C block:\n");

       //MatView(Schur->Compl->A,PETSC_VIEWER_STDOUT_WORLD);
       MatAXPY(Compl->get_matrix(), 1, xA, SUBSET_NONZERO_PATTERN);
       DBGMSG("C block:\n");
       //MatView(Schur->Compl->A,PETSC_VIEWER_STDOUT_WORLD);

    }

    form_rhs();
    state=formed;

}

void SchurComplement::form_rhs()
{
    PetscScalar *rhs_interior, *rhs_boundary, *rhs_update_array;
    PetscInt *ind_interior;
    PetscInt *interface_subdomain_indices;
    PetscInt size, loc_size;
    Vec rhs1_vec, rhs1_multiplied_vec;
    Vec RHS2_update;
    Vec rhs_update_big;
    ISLocalToGlobalMapping B_map;
    int i;
    int ind, ind_global;
    int orig_lsize, locSizeB_vec;

    // compute the SchurRHS
    if      (Orig->type == LinSys::MAT_IS)
    {
       VecGetArray(RHS1,&rhs_interior);

       //for (i = 0; i<locSizeA; i++) 
       //{
       //   xprintf(Msg,"prvek %d : %f\n",i,rhs_interior[i]);
       //}

       VecCreateSeqWithArray(PETSC_COMM_SELF,locSizeA,rhs_interior,&rhs1_vec);
       VecRestoreArray(RHS1,PETSC_NULL);

       VecCreateSeq(PETSC_COMM_SELF,locSizeB,&rhs1_multiplied_vec);
       MatMultTranspose(IAB_sub,rhs1_vec,rhs1_multiplied_vec);

       VecGetArray(rhs1_multiplied_vec,&rhs_boundary);
       //VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB,PETSC_DETERMINE,rhs_boundary,&RHS2_update);

       VecCreateMPI(PETSC_COMM_WORLD,Orig->vec_lsize(),PETSC_DETERMINE,&rhs_update_big);
       //for (i = 0; i<locSizeB; i++) 
       //{
       //   xprintf(Msg,"prvek %d : %f\n",i,rhs_boundary[i]);
       //}

       for (i = 0; i<locSizeB; i++)
       {
	  // local index on subdomain
	  ind = locSizeA + i;

	  // global index
	  ind_global = Orig->subdomain_indices[ind];

          VecSetValue(rhs_update_big,ind_global,rhs_boundary[i],INSERT_VALUES);

       }
       VecAssemblyBegin(rhs_update_big);
       VecAssemblyEnd(rhs_update_big);

       // create RHS sub vecs RHS1, RHS2
       orig_lsize = Orig->vec_lsize();
       locSizeB_vec = orig_lsize - locSizeA;
       VecGetArray(rhs_update_big,&rhs_update_array);
       VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB_vec,PETSC_DETERMINE,rhs_update_array+locSizeA,&(RHS2_update));
       VecRestoreArray(rhs_update_big,PETSC_NULL);

       VecAXPY(Compl->get_rhs(),1,RHS2_update);
       VecAXPY(Compl->get_rhs(),-1,RHS2);

       VecDestroy(rhs1_vec);
       VecDestroy(rhs1_multiplied_vec);
       VecDestroy(rhs_update_big);
       VecDestroy(RHS2_update);
    }
    else if (Orig->type == LinSys::MAT_MPIAIJ)
    {
       MatMultTranspose(IAB,RHS1,Compl->get_rhs());
       VecAXPY(Compl->get_rhs(),-1,RHS2);
    }

}

/**
 *  @brief Scale formed complement system. Mainly to make it positive definite.
 */

void SchurComplement ::scale(double scalar)
{
    ASSERT( state == formed, "Object in wrong state!\n");
    MatScale(Compl->get_matrix(), scalar);
    VecScale(Compl->get_rhs(), scalar);
}

/**
 * COMPUTE ELIMINATED PART OF THE ORIG. SYS. & RESTORE RHS and SOLUTION VECTORS
 *  x_1 = IA * RHS_1 - IAB * x_2
 */

void SchurComplement::resolve()
{
    F_ENTRY;

    MatMult(IAB,Compl->get_solution(),Sol1);
    VecScale(Sol1,-1);
    MatMultAdd(IA,RHS1,Sol1,Sol1);
}

/**
 * SCHUR COMPLEMENT destructor
 */
SchurComplement :: ~SchurComplement() {

    F_ENTRY;

    MatDestroy(IA);
    MatDestroy(IA_sub);
    MatDestroy(B);
    MatDestroy(B_sub);
    MatDestroy(Bt);
    MatDestroy(Bt_sub);
    MatDestroy(xA);
    MatDestroy(xA_sub);
    MatDestroy(IAB);
    MatDestroy(IAB_sub);
    ISDestroy(IsA);
    ISDestroy(IsB);
    ISDestroy(fullIsA);
    ISDestroy(fullIsB);
    VecDestroy(RHS1);
    VecDestroy(RHS2);
    VecDestroy(Sol1);
    VecDestroy(Sol2);
    delete Compl;
}
