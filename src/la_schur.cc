/*!
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
 * TODO:
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
    int i,n;
    IS ISBloc;

    F_ENTRY;

    // get distribution of original marix
    int orig_first, orig_lsize;
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

    // indicate first construction
    mat_reuse=MAT_INITIAL_MATRIX;;

    DBGMSG("ISs :\n");
    //ISView(Schur->IsA,PETSC_VIEWER_STDOUT_WORLD);
    //ISView(Schur->IsB,PETSC_VIEWER_STDOUT_WORLD);
    //ISView(Schur->fullIsA,PETSC_VIEWER_STDOUT_SELF);
    //ISView(Schur->fullIsB,PETSC_VIEWER_STDOUT_SELF);


    // create complement system
    // TODO: introduce LS as true object, clarify its internal states
    // create RHS sub vecs RHS1, RHS2
    VecGetArray(Orig->get_rhs(),&rhs_array);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeA,PETSC_DETERMINE,rhs_array,&(RHS1));
    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB,PETSC_DETERMINE,rhs_array+locSizeA,&(RHS2));
    VecRestoreArray(Orig->get_rhs(),PETSC_NULL);

    // create Solution sub vecs Sol1, Compl->solution
    VecGetArray(Orig->get_solution(),&sol_array);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeA,PETSC_DETERMINE,sol_array,&(Sol1));
    Compl = new LinSys_MPIAIJ(locSizeB, sol_array+locSizeA);

    VecCreateMPIWithArray(PETSC_COMM_WORLD,locSizeB,PETSC_DETERMINE,sol_array+locSizeA,&(Sol2));
    VecRestoreArray(Orig->get_solution(),PETSC_NULL);

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
    ASSERT(state == created, "Object in wrong state.\n");

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
    MatMatMult(IA, B, mat_reuse, 1.1 ,&(IAB)); // 6/7 - fill estimate
    DBGMSG(" IAB:\n");
    //MatView(Schur->IAB,PETSC_VIEWER_STDOUT_WORLD);
    // compute xA=Bt* IAB = Bt * IA * B
    MatGetSubMatrix(Orig->get_matrix(), IsB, fullIsA, locSizeA, mat_reuse, &(Bt));
    MatMatMult(Bt, IAB, mat_reuse, 1.5 ,&(xA)); // 1.1 - fill estimate
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

    // compute the SchurRHS
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
    MatDestroy(B);
    MatDestroy(Bt);
    MatDestroy(xA);
    MatDestroy(IAB);
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
