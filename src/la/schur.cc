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
 * @ingroup la
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
#include <algorithm>
#include <limits>
#include <petscmat.h>
#include <armadillo>

#include "la/distribution.hh"
#include "la/local_to_global_map.hh"
#include "solve.h"
#include "system/system.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"

/**
 *  Create Schur complement system.
 *  @param[in] orig  : original system
 *  @param[in] inv_a : inversion of the A block
 *  @param[in] ia    : index set of the A block,
 *   default continuous given by inv_a:
 *   proc    1      2      3
 *
 *   Orig:   ****** ****** ****
 *   IA  :   ***    **     ***
 *
  *
 */

SchurComplement :: SchurComplement(LinSys *orig, Mat & inv_a, IS ia)
: IA(inv_a), IsA(ia), state(created), Orig(orig)

{
    PetscScalar *rhs_array, *sol_array;
    //const PetscInt *IsAIndices;
    //PetscErrorCode ierr;

    //int i;

    int orig_first;

    // MATIS vars
    //Mat orig_mat_sub;
    //VecScatter ScatterToA;
    //VecScatter ScatterToB;

    //int locSizeB_vec;
    //int vec_orig_low,vec_orig_high;

    // initialize variables
    IA_sub  = NULL;
    B       = NULL; 
    Bt      = NULL;        
    B_sub   = NULL;
    Bt_sub  = NULL; 
    xA      = NULL; 
    xA_sub  = NULL;
    IAB     = NULL; 
    IAB_sub = NULL;        
    sub_vec_block2 = NULL;        

    F_ENTRY;

    /* // check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {
       // it is assumed that Schur complement may be formed block-wise

       // get local submatrix of A_inverse
       ierr = MatISGetLocalMat(IA, &IA_sub);
       ASSERT(ierr == 0,"Error in MatISGetLocalMat.");

       // get local submatrix of original system
       ierr = MatISGetLocalMat(Orig->get_matrix(), &orig_mat_sub);
       ASSERT(ierr == 0,"Error in MatISGetLocalMat.");

       // find original size
       PetscInt m,n;
       ierr = MatGetSize(orig_mat_sub, &m, &n);
       ASSERT(m == n,"Assumed square matrix.");

       // set size
       orig_sub_size = m;
       //DBGMSG("orig block size %d",orig_sub_size);

       // find inverted block size
       ierr = MatGetSize(IA_sub, &m, &n);
       ASSERT(m == n,"Assumed square matrix.");

       // set size
       loc_size_A = m;
       //DBGMSG("A block size %d",locSizeA);

       // size of B block
       locSizeB = orig_sub_size - loc_size_A;

       VecGetOwnershipRange(Orig->get_rhs(),&vec_orig_low,&vec_orig_high);

       if (IsA == NULL) {
           // create A block index set
           // assume 'a_inv->local_size' be local part of A block
           ISCreateStride(PETSC_COMM_SELF,loc_size_A,vec_orig_low,1,&IsA);

       } 
       // obtain index set local to subdomain
       ISGetIndices(IsA,&IsAIndices);
       IsALocalIndices = new PetscInt[loc_size_A];
       int shift = IsAIndices[0];
       for (i = 0; i < loc_size_A; i++) {
	   IsALocalIndices[i] = IsAIndices[i] - shift;
       }
       ISRestoreIndices(IsA,&IsAIndices);

       ISCreateGeneral(PETSC_COMM_SELF,loc_size_A,IsALocalIndices,PETSC_USE_POINTER,&IsA_sub);
       //DBGPRINT_INT("pole_lokalnich_indexu",locSizeA,IsALocalIndices);
       MPI_Barrier(PETSC_COMM_WORLD);

       ISAllGather(IsA,&fullIsA);

       // create B block index set
       orig_lsize = Orig->vec_lsize();
       locSizeB_vec = orig_lsize - loc_size_A;

       ISComplement(IsA, vec_orig_low, vec_orig_high, &IsB);
       ISAllGather(IsB,&fullIsB);

       // index set local to subdomain
       ISComplement(IsA_sub, 0, orig_sub_size, &IsB_sub);


       //DBGMSG("ISs :\n");
       //ISView(IsA,PETSC_VIEWER_STDOUT_WORLD);
       //ISView(IsB,PETSC_VIEWER_STDOUT_WORLD);
       //ISView(fullIsA,PETSC_VIEWER_STDOUT_SELF);
       //ISView(fullIsB,PETSC_VIEWER_STDOUT_SELF);

       // create complement system
       // TODO: introduce LS as true object, clarify its internal states
       // create RHS sub vecs RHS1, RHS2
       ierr = VecCreateMPI(PETSC_COMM_WORLD,loc_size_A,PETSC_DETERMINE,&(RHS1));
       ierr = VecCreateMPI(PETSC_COMM_WORLD,locSizeB_vec,PETSC_DETERMINE,&(RHS2));

       // create scatters
       ierr = VecScatterCreate(Orig->get_rhs( ), IsA, RHS1, PETSC_NULL, &ScatterToA);
       ierr = VecScatterCreate(Orig->get_rhs( ), IsB, RHS2, PETSC_NULL, &ScatterToB);

       // create Solution sub vecs Sol1, Compl->solution
       ierr = VecCreateMPI(PETSC_COMM_WORLD, loc_size_A,     PETSC_DETERMINE, &(Sol1));
       ierr = VecCreateMPI(PETSC_COMM_WORLD, locSizeB_vec, PETSC_DETERMINE, &(Sol2));

       // Prepare local data
       ierr = VecScatterBegin( ScatterToA, Orig->get_rhs(),      RHS1, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterBegin( ScatterToB, Orig->get_rhs(),      RHS2, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterEnd  ( ScatterToA, Orig->get_rhs(),      RHS1, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterEnd  ( ScatterToB, Orig->get_rhs(),      RHS2, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterBegin( ScatterToA, Orig->get_solution(), Sol1, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterBegin( ScatterToB, Orig->get_solution(), Sol2, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterEnd  ( ScatterToA, Orig->get_solution(), Sol1, INSERT_VALUES, SCATTER_FORWARD);
       ierr = VecScatterEnd  ( ScatterToB, Orig->get_solution(), Sol2, INSERT_VALUES, SCATTER_FORWARD);

       // Destroy Scatters
       ierr = VecScatterDestroy ( &(ScatterToA) );
       ierr = VecScatterDestroy ( &(ScatterToB) );

    } */

    // check type of LinSys
    if (typeid(*Orig) == typeid(LinSys_PETSC)){
       // get distribution of original marix
       MatGetOwnershipRange(Orig->get_matrix(),&orig_first,PETSC_NULL);
       MatGetLocalSize(Orig->get_matrix(),&orig_lsize,PETSC_NULL);

       // create A block index set
       if (IsA == NULL) {
           // assume 'a_inv->local_size' be local part of A block
           MatGetLocalSize(IA,&loc_size_A,PETSC_NULL);
           ISCreateStride(PETSC_COMM_WORLD,loc_size_A,orig_first,1,&IsA);
       } else {
           ISGetLocalSize(IsA, &loc_size_A);
       }
       ISAllGather(IsA,&fullIsA);

       // create B block index set
       locSizeB = orig_lsize-loc_size_A;
       ISCreateStride(PETSC_COMM_WORLD,locSizeB,orig_first+loc_size_A,1,&IsB);
       ISAllGather(IsB,&fullIsB);


       //DBGMSG("ISs :\n");
       //ISView(IsA,PETSC_VIEWER_STDOUT_WORLD);
       //ISView(IsB,PETSC_VIEWER_STDOUT_WORLD);
       //ISView(fullIsA,PETSC_VIEWER_STDOUT_SELF);
       //ISView(fullIsB,PETSC_VIEWER_STDOUT_SELF);


       // create complement system
       // TODO: introduce LS as true object, clarify its internal states
       // create RHS sub vecs RHS1, RHS2
       VecGetArray(Orig->get_rhs(),&rhs_array);
       VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,rhs_array,&(RHS1));

       // create Solution sub vecs Sol1, Compl->solution
       VecGetArray(Orig->get_solution(),&sol_array);
       VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,sol_array,&(Sol1));

       VecCreateMPIWithArray(PETSC_COMM_WORLD,1,locSizeB,PETSC_DETERMINE,rhs_array+loc_size_A,&(RHS2));
       VecCreateMPIWithArray(PETSC_COMM_WORLD,1,locSizeB,PETSC_DETERMINE,sol_array+loc_size_A,&(Sol2));

       VecRestoreArray(Orig->get_rhs(),&rhs_array);
       VecRestoreArray(Orig->get_solution(),&sol_array);

       VecGetArray( Sol2, &sol_array );
       ds_ = new Distribution(locSizeB, PETSC_COMM_WORLD);
       Compl = new LinSys_PETSC( ds_, PETSC_COMM_WORLD );
       Compl->set_solution(sol_array);
       VecRestoreArray( Sol2, &sol_array );

    } else if (typeid(*Orig) == typeid(LinSys_BDDC)){
    	xprintf(Warn, "Schur complement for LinSys::BDDC type shouldn't be constructed!\n");
    }

    /* // check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {

       const PetscInt *rangesAblock;
       VecGetOwnershipRanges(RHS1,&rangesAblock);

       Distribution new_ds(locSizeB, PETSC_COMM_WORLD);
       boost::shared_ptr<LocalToGlobalMap> global_row_4_sub_row_new;
       global_row_4_sub_row_new=boost::make_shared<LocalToGlobalMap>(new_ds);

       int indB = 0;
       int subInd, indglb;
       int myid;
       int proc;
       int shift;
       int new_index;
       // pick indices outside interior block
       myid = Orig->vec_ds.myp( );

       const PetscInt *IsBLocalIndices;

       ISGetIndices(IsB_sub, &IsBLocalIndices);

       for (int i = 0; i < locSizeB; i++)
       { 

           subInd = IsBLocalIndices[i];

          indglb = Orig->subdomain_indices[subInd];

          // find processor of this value
          proc = Orig->vec_ds.get_proc(indglb);

          //get shifted number
          shift = rangesAblock[proc+1]-1;
          new_index = indglb - shift - 1;

          global_row_4_sub_row_new->insert(new_index); // vec_orig_last is returned
	                                                  // larger of one by PETSc
          indB = indB + 1;
       }

       ISRestoreIndices(IsB_sub, &IsBLocalIndices);

       if (indB != locSizeB) 
          xprintf(Err,"Length of second block in Schur complement mismatch: %d,%d \n",indB, locSizeB);
	  

       //DBGPRINT_INT("pole_indexu_nove",locSizeB,global_row_4_sub_row_new);

       VecGetArray( Sol2, &sol_array );

       Compl = new LinSys_MATIS( global_row_4_sub_row_new, &(sol_array[0]) );
       Compl->start_allocation();

       VecRestoreArray( Sol2, &sol_array );

    }*/


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
}

SchurComplement :: SchurComplement(LinSys *orig, IS ia, PetscInt max_size_submat)
: IsA(ia), state(created), Orig(orig)
{
        xprintf(Msg, "Constructor SchurComplement\n");

        PetscInt m, n;
        PetscErrorCode ierr;
        PetscScalar *rhs_array, *sol_array;
        int orig_first;

        int n_proc, rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &n_proc);

        // check type of LinSys, dimensions of matrix, index set
        ierr = MatGetSize(Orig->get_matrix(), &m, &n);
        ASSERT(typeid(*Orig) == typeid(LinSys_PETSC), "Assumed PETSC type of Orig object.\n");
        //ASSERT(Orig->type == LinSys::PETSC, "Assumed PETSC type of Orig object.\n");
        ASSERT(m == n, "Assumed square matrix.\n" );
        ASSERT(IsA != NULL, "Index set IsA is not defined.\n" );

        // initialize variables
        IA_sub  = NULL;
        B       = NULL;
        Bt      = NULL;
        B_sub   = NULL;
        Bt_sub  = NULL;
        xA      = NULL;
        xA_sub  = NULL;
        IAB     = NULL;
        IAB_sub = NULL;
        sub_vec_block2 = NULL;

        F_ENTRY;

        {
           // get distribution of original matrix
           MatGetOwnershipRange(Orig->get_matrix(),&orig_first,PETSC_NULL);
           MatGetLocalSize(Orig->get_matrix(),&orig_lsize,PETSC_NULL);
           MatGetSubMatrix(Orig->get_matrix(), IsA, IsA, MAT_INITIAL_MATRIX, &IA);
           DBGMSG(" - rank: %d, size: %d, orig_first: %d, orig_lsize: %d \n", rank, m, orig_first, orig_lsize);
           //MatView(IA,PETSC_VIEWER_STDOUT_WORLD);

           // create A block index set
           ISGetLocalSize(IsA, &loc_size_A);
           ISAllGather(IsA,&fullIsA);
           DBGMSG(" - rank: %d, locSizeA: %d\n", rank, loc_size_A);
           //ISView(IsA, PETSC_VIEWER_STDOUT_WORLD);

           // create B block index set
           locSizeB = orig_lsize-loc_size_A;
           ISCreateStride(PETSC_COMM_WORLD,locSizeB,orig_first+loc_size_A,1,&IsB);
           ISAllGather(IsB,&fullIsB);
           DBGMSG(" - rank: %d, locSizeB: %d\n", rank, locSizeB);
           //ISView(IsB, PETSC_VIEWER_STDOUT_WORLD);

           // create complement system
           // TODO: introduce LS as true object, clarify its internal states
           // create RHS sub vecs RHS1, RHS2
           VecGetArray(Orig->get_rhs(),&rhs_array);
           VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,rhs_array,&(RHS1));

           // create Solution sub vecs Sol1, Compl->solution
           VecGetArray(Orig->get_solution(),&sol_array);
           VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,sol_array,&(Sol1));

           VecCreateMPIWithArray(PETSC_COMM_WORLD,1,locSizeB,PETSC_DETERMINE,rhs_array+loc_size_A,&(RHS2));
           VecCreateMPIWithArray(PETSC_COMM_WORLD,1,locSizeB,PETSC_DETERMINE,sol_array+loc_size_A,&(Sol2));

           VecRestoreArray(Orig->get_rhs(),&rhs_array);
           VecRestoreArray(Orig->get_solution(),&sol_array);

           VecGetArray( Sol2, &sol_array );
           ds_ = new Distribution(locSizeB, PETSC_COMM_WORLD);
           Compl = new LinSys_PETSC( ds_, PETSC_COMM_WORLD );
           Compl->set_solution(sol_array);
           VecRestoreArray( Sol2, &sol_array );
        }

    	PetscInt ncols, pos_start, pos_start_IA;
    	MatGetOwnershipRange(Orig->get_matrix(),&pos_start,PETSC_NULL);
    	MatGetOwnershipRange(IA,&pos_start_IA,PETSC_NULL);
        DBGMSG(" - rank: %d, pos_start_IA: %d, pos_start: %d, loc_size_A: %d \n", rank, pos_start_IA, pos_start, loc_size_A);

        std::vector<PetscInt> submat_rows;
        const PetscInt *cols;
        const PetscScalar *vals;

        std::vector<unsigned int> processed_rows(loc_size_A,0);

        unsigned int mat_block=1;   //actual processed block of matrix
        for(unsigned int loc_row=0; loc_row < processed_rows.size(); loc_row++) {
            if (processed_rows[loc_row] != 0) continue;

            	PetscInt min=std::numeric_limits<int>::max(), max=-1, size_submat;
            	unsigned int b_vals = 0; // count of values stored in B-block of Orig system
                submat_rows.clear();
                ierr = MatGetRow(Orig->get_matrix(), loc_row + pos_start, &ncols, &cols, PETSC_NULL);
                for (PetscInt i=0; i<ncols; i++) {
                	if (cols[i] >= pos_start+loc_size_A) {
                		b_vals++;
                	} else if (cols[i] < min) {
                		min=cols[i];
                	} else if (cols[i] > max) {
                		max=cols[i];
                	}
                }
                size_submat = max - min + 1;
                ASSERT(ncols-b_vals == size_submat, "Submatrix cannot contains empty values.\n");
                ASSERT_LESS(size_submat , max_size_submat);

                ierr = MatRestoreRow(Orig->get_matrix(), loc_row + pos_start, &ncols, &cols, PETSC_NULL);
                arma::mat submat2(size_submat, size_submat);
                submat2.zeros();
                for (PetscInt i=0; i<size_submat; i++) {
                	processed_rows[ loc_row + i ] = mat_block;
                	submat_rows.push_back( i + loc_row + pos_start_IA );
                    ierr = MatGetRow(Orig->get_matrix(), i + loc_row + pos_start, &ncols, &cols, &vals);
                    for (PetscInt j=0; j<ncols; j++) {
                    	if (cols[j] < pos_start+loc_size_A) submat2( i, cols[j] - loc_row - pos_start ) = vals[j];
                    }
                    ierr = MatRestoreRow(Orig->get_matrix(), i + loc_row + pos_start, &ncols, &cols, &vals);
                }
                // test output
//                xprintf(Msg, "__ Get submat: rank %d, MIN-MAX %d %d, size %d\n", rank, min, max, size_submat);
//                for (int i=0; i<size_submat; i++) {
//                    for (int j=0; j<size_submat; j++) {
//                        xprintf(Msg, "%2.0f ", submat2(i,j));
//                    }
//                    xprintf(Msg, "\n");
//                }
//                xprintf(Msg, "\n");

                // get inversion matrix
                arma::mat invmat = submat2.i();
                // stored to inversion IA matrix
                const PetscInt* rows = &submat_rows[0];
                MatSetValues(IA, submat_rows.size(), rows, submat_rows.size(), rows, invmat.memptr(), INSERT_VALUES);

                mat_block++;
        }

        MatAssemblyBegin(IA, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(IA, MAT_FINAL_ASSEMBLY);
}


#if 0
/**
 * @brief Form Schur complement. Call solve. Resolve original solution.
 * TODO: should be better when LinSys is full object.
 *
 */

void SchurComplement::solve(Solver *solver)
{
    if (state != formed) form_schur();
    ASSERT(state == formed, "Object in wrong state.\n" );

    solve_system(solver, get_system() );

    resolve();
    state=solved;

}
#endif

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
    PetscErrorCode ierr = 0;
    //Mat orig_mat_sub;
    //Mat local_compl_aux;
    MatReuse mat_reuse;        // reuse structures after first computation of schur

    mat_reuse=MAT_REUSE_MATRIX;
    if (state==created) mat_reuse=MAT_INITIAL_MATRIX; // indicate first construction

    /*// check type of LinSys
    if      (Orig->type == LinSys::MAT_IS)
    {
       // get local submatrix of original system
       ierr = MatISGetLocalMat(Orig->get_matrix(), &orig_mat_sub);
       ASSERT(ierr == 0,"Error in MatISGetLocalMat.");

       // B, locSizeB removed
       ierr+=MatGetSubMatrix(orig_mat_sub, IsA_sub, IsB_sub, mat_reuse, &B_sub);

       // A^-1 * B
       ierr+=MatMatMult(IA_sub, B_sub, mat_reuse, 1.0 ,&(IAB_sub)); // 6/7 - fill estimate

       // B^T,  locSizeA removed
       ierr+=MatGetSubMatrix(orig_mat_sub, IsB_sub, IsA_sub, mat_reuse, &Bt_sub);

       // B^T*A^-1*B
       ierr+=MatMatMult(Bt_sub, IAB_sub, mat_reuse, 1.9 ,&(xA_sub)); // 1.1 - fill estimate (PETSC report values over 1.8)

       // get C block TODO: matrix reuse,  locSizeB removed
       ierr+=MatGetSubMatrix(orig_mat_sub, IsB_sub, IsB_sub, MAT_INITIAL_MATRIX, &local_compl_aux);
       ierr+=MatCopy(local_compl_aux, Compl->local_matrix,DIFFERENT_NONZERO_PATTERN);
       ierr+=MatDestroy(&local_compl_aux);

       // compute complement = (-1)cA+xA = Bt*IA*B - C
       ierr+=MatScale(Compl->local_matrix,-1.0);
       ierr+=MatAXPY(Compl->local_matrix, 1, xA_sub, DIFFERENT_NONZERO_PATTERN);

       // assemble final MATIS matrix
       ierr+=MatAssemblyBegin(Compl->get_matrix(),MAT_FINAL_ASSEMBLY);
       ierr+=MatAssemblyEnd(Compl->get_matrix(),MAT_FINAL_ASSEMBLY);

       ASSERT( ierr == 0, "PETSC Error during claculation of Schur complement.\n");

    }
    else if      (Orig->type == LinSys::MAT_MPIAIJ)*/
    {

       //DBGMSG("Compute Schur complement of\n");
       //MatView(Orig->get_matrix(),PETSC_VIEWER_STDOUT_WORLD);
       //DBGMSG("inverse IA:\n");
       //MatView(Schur->IA,PETSC_VIEWER_STDOUT_WORLD);
       // compose Schur complement
       // Petsc need some fill estimate for results of multiplication in form nnz(A*B)/(nnz(A)+nnz(B))
       // for the first Schur compl: IA*B is bounded by ( d*(d+1) )/( d*d+2*d ) <= 5/6 for d<=4
       //                            B'*IA*B bounded by ( (d+1)*(d+1) )/ ( d*(d+1) + d ) ~ 1
       // for the second Schur :      IA*B have fill ratio ~ 1.
       //                            B'*IA*B  ...         ( N/2 *(2*N-1) )/( 2 + 2*N ) <= 1.4
       // nevertheless Petsc does not allows fill ratio below 1. so we use 1.1 for the first
       // and 1.5 for the second multiplication

       // compute IAB=IA*B, locSizeB removed
       ierr+=MatGetSubMatrix(Orig->get_matrix(), IsA, IsB, mat_reuse, &B);
       //DBGMSG(" B:\n");
       //MatView(Schur->B,PETSC_VIEWER_STDOUT_WORLD);
       ierr+=MatMatMult(IA, B, mat_reuse, 1.0 ,&(IAB)); // 6/7 - fill estimate
       //DBGMSG(" IAB:\n");
       //MatView(Schur->IAB,PETSC_VIEWER_STDOUT_WORLD);
       // compute xA=Bt* IAB = Bt * IA * B, locSizeA removed
       ierr+=MatGetSubMatrix(Orig->get_matrix(), IsB, IsA, mat_reuse, &(Bt));
       ierr+=MatMatMult(Bt, IAB, mat_reuse, 1.9 ,&(xA)); // 1.1 - fill estimate (PETSC report values over 1.8)
       //DBGMSG("xA:\n");
       //MatView(Schur->xA,PETSC_VIEWER_STDOUT_WORLD);

       // get C block, locSizeB removed
       ierr+=MatGetSubMatrix( Orig->get_matrix(), IsB, IsB, mat_reuse, const_cast<Mat *>( &(Compl->get_matrix()) ) );
       // compute complement = (-1)cA+xA = Bt*IA*B - C
       ierr+=MatScale(Compl->get_matrix(),-1.0);
       //DBGMSG("C block:\n");

       //MatView(Schur->Compl->A,PETSC_VIEWER_STDOUT_WORLD);
       ierr+=MatAXPY(Compl->get_matrix(), 1, xA, SUBSET_NONZERO_PATTERN);
       //DBGMSG("C block:\n");
       //MatView(Schur->Compl->A,PETSC_VIEWER_STDOUT_WORLD);
       //
       //TODO: MatAXPY - umoznuje nasobit -1, t.j. bylo by lepe vytvorit konvencni Schuruv doplnek zde,
       // a ve funkci schur1 (a ne v schur2) uzit metodu "scale" z tohoto objektu - kvuli negativni definitnosti
       // usetri se tim jeden MatScale

       ASSERT( ierr == 0, "PETSC Error during calculation of Schur complement.\n");
    }

    /*
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matAinv.output",&myViewer);
    PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
    MatView( IA, myViewer );
    PetscViewerDestroy(myViewer);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matSchur.output",&myViewer);
    PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
    MatView( Compl->get_matrix( ), myViewer );
    PetscViewerDestroy(myViewer);
*/

    form_rhs();

    state=formed;

}

void SchurComplement::form_rhs()
{
    //PetscScalar *rhs_interior, *subdomain_rhs_array;
    //PetscErrorCode ierr;
    //const PetscInt *IsALocalIndices;
    //Vec rhs1_vec;
    //Vec RHS2_update;
    //Vec subdomain_rhs;
    //int i;
    //int orig_lsize, locSizeB_vec;
    //int locIndex;

    // compute the SchurRHS
    /*if      (Orig->type == LinSys::MAT_IS)
    {

       VecCreateSeq(PETSC_COMM_SELF, orig_sub_size, &subdomain_rhs);

       LinSys_MATIS *ls_IS_Orig = dynamic_cast<LinSys_MATIS*>(Orig);
       ierr = VecScatterBegin(ls_IS_Orig->get_scatter(),Orig->get_rhs(),subdomain_rhs,INSERT_VALUES,SCATTER_FORWARD);
       ierr = VecScatterEnd(ls_IS_Orig->get_scatter(),Orig->get_rhs(),subdomain_rhs,INSERT_VALUES,SCATTER_FORWARD);

       ISGetIndices( IsA_sub, &IsALocalIndices );
       VecGetArray( subdomain_rhs, &subdomain_rhs_array );

       rhs_interior = new PetscScalar[loc_size_A];

       for (i = 0; i<loc_size_A; i++)
       {
	  locIndex = IsALocalIndices[i];
          rhs_interior[i] = subdomain_rhs_array[locIndex];
       }
       VecRestoreArray( subdomain_rhs, &subdomain_rhs_array );
       ISRestoreIndices( IsA_sub, &IsALocalIndices );


       VecCreateSeqWithArray(PETSC_COMM_SELF,1, loc_size_A,rhs_interior,&rhs1_vec);

       VecCreateSeq(PETSC_COMM_SELF,locSizeB,&sub_vec_block2);
       MatMultTranspose(IAB_sub,rhs1_vec,sub_vec_block2);

       // create RHS sub vecs RHS1, RHS2
       orig_lsize   = Orig->vec_lsize();
       locSizeB_vec = orig_lsize - loc_size_A;
       VecCreateMPI(PETSC_COMM_WORLD, locSizeB_vec,PETSC_DETERMINE,&(RHS2_update));

       LinSys_MATIS *ls_IS_Compl = dynamic_cast<LinSys_MATIS*>(Compl);
       ierr = VecScatterBegin( ls_IS_Compl->get_scatter(), sub_vec_block2, RHS2_update,  ADD_VALUES, SCATTER_REVERSE);
       ierr = VecScatterEnd(   ls_IS_Compl->get_scatter(), sub_vec_block2, RHS2_update,  ADD_VALUES, SCATTER_REVERSE);

       // Prepare reduced RHS
       VecSet(Compl->get_rhs( ),0.0);
       VecAXPY(Compl->get_rhs( ),1,RHS2_update);
       VecAXPY(Compl->get_rhs( ),-1,RHS2);

       VecDestroy(&(rhs1_vec));
       delete [ ] rhs_interior;
       VecDestroy(&(RHS2_update));

    else if (Orig->type == LinSys::MAT_MPIAIJ)*/

    MatMultTranspose(IAB,RHS1,Compl->get_rhs());
    VecAXPY(Compl->get_rhs(),-1,RHS2);

    state=formed;
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

    //PetscErrorCode ierr;
    //Vec sol1_vec_loc;
    //PetscScalar *sol1_array_loc;

    /*if      (Compl->type == LinSys::MAT_IS)
    {
       //  scatter the global vector x into the local work vector
       LinSys_MATIS *ls_IS = dynamic_cast<LinSys_MATIS*>(Compl);
       ierr = VecScatterBegin(ls_IS->get_scatter(),Compl->get_solution(),sub_vec_block2, INSERT_VALUES,SCATTER_FORWARD);
       ierr = VecScatterEnd(ls_IS->get_scatter(),Compl->get_solution(),  sub_vec_block2, INSERT_VALUES,SCATTER_FORWARD);

       VecGetArray(Sol1,&sol1_array_loc);
       VecCreateSeqWithArray(PETSC_COMM_SELF,1, loc_size_A,sol1_array_loc,&sol1_vec_loc);

       MatMult(IAB_sub,sub_vec_block2,sol1_vec_loc);
       VecRestoreArray(Sol1,&sol1_array_loc);

       VecDestroy(&(sol1_vec_loc));

    }
    else if (Orig->type == LinSys::MAT_MPIAIJ) // */
    {
       MatMult(IAB,Compl->get_solution(),Sol1);
    }

    VecScale(Sol1,-1);

    MatMultAdd(IA,RHS1,Sol1,Sol1);

    // join local portions of solution 
    /*if      (Orig->type == LinSys::MAT_IS)
    {
        PetscScalar * sol_array_loc;
        VecGetArray( Orig->get_solution(),&sol_array_loc );

        PetscScalar * sol2_array_loc;
        VecGetArray( Sol1, &sol1_array_loc );
        VecGetArray( Compl->get_solution(), &sol2_array_loc );

        for (int i = 0; i<loc_size_A; i++)
        {
            sol_array_loc[i] = sol1_array_loc[i];
        }
        int locSizeB_vec = orig_lsize - loc_size_A;
        for (int i = 0; i<locSizeB_vec; i++) 
        {
            int index = loc_size_A + i;
            sol_array_loc[ index ] = sol2_array_loc[i];
        }

        VecRestoreArray( Orig->get_solution(),  &sol_array_loc );
        VecRestoreArray( Sol1, &sol1_array_loc );
        VecRestoreArray( Compl->get_solution(), &sol2_array_loc );
    } // */

}

/**
 * SCHUR COMPLEMENT destructor
 */
SchurComplement :: ~SchurComplement() {

    F_ENTRY;

    if ( B  != NULL )             MatDestroy(&B);
    if ( Bt != NULL )             MatDestroy(&Bt);
    if ( xA != NULL )             MatDestroy(&xA);
    if ( IAB != NULL )            MatDestroy(&IAB);
    if ( IsA != NULL )            ISDestroy(&IsA);
    if ( IsB != NULL )            ISDestroy(&IsB);
    if ( fullIsA != NULL )        ISDestroy(&fullIsA);
    if ( fullIsB != NULL )        ISDestroy(&fullIsB);
    if ( RHS1 != NULL )           VecDestroy(&RHS1);
    if ( RHS2 != NULL )           VecDestroy(&RHS2);
    if ( Sol1 != NULL )           VecDestroy(&Sol1);
    if ( Sol2 != NULL )           VecDestroy(&Sol2);
    if ( B_sub != NULL )          MatDestroy(&B_sub);
    if ( Bt_sub != NULL )         MatDestroy(&Bt_sub);
    if ( xA_sub != NULL )         MatDestroy(&xA_sub);
    if ( IAB_sub != NULL )        MatDestroy(&IAB_sub);
    if ( sub_vec_block2 != NULL ) VecDestroy(&sub_vec_block2);

//    if      (Orig->type == LinSys::MAT_IS) {
//        delete [] IsALocalIndices;
//    }
    delete Compl;
}
