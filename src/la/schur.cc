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
 * @file    schur.cc
 * @ingroup la
 * @brief   Assembly explicit Schur complement for the given linear system.
 *          Provides method for resolution of the full original vector of unknowns.
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
#include <petscis.h>

#include "system/sys_profiler.hh"
#include "la/distribution.hh"
#include "la/local_to_global_map.hh"
#include "system/system.hh"
#include "la/linsys.hh"
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

SchurComplement::SchurComplement(Distribution *ds, IS ia, IS ib)
: LinSys_PETSC(ds), IsA(ia), IsB(ib), state(created)
{
        // check index set
        OLD_ASSERT(IsA != NULL, "Index set IsA is not defined.\n" );

        // initialize variables
        Compl   = NULL;
        A       = NULL;
        IA      = NULL;
        B       = NULL;
        Bt      = NULL;
        C       = NULL;
        xA      = NULL;
        IAB     = NULL;
        IsB     = NULL;
        RHS1    = NULL;
        RHS2    = NULL;
        Sol1    = NULL;
        Sol2    = NULL;
        rhs1sc  = NULL;
        rhs2sc  = NULL;
        sol1sc  = NULL;
        sol2sc  = NULL;
        ds_     = NULL;

        // create A block index set
        ISGetLocalSize(IsA, &loc_size_A);

        // create B block index set
        if ( IsB == nullptr ) {

            // compute dual indexset

            vector<bool> a_used(ds->lsize(), false);

            const PetscInt *loc_indices;
            ISGetIndices(IsA, &loc_indices);
            for(int i=0; i < loc_size_A; i++)
                a_used[ loc_indices[i] - ds->begin()] = true;
            ISRestoreIndices(IsA, &loc_indices);

            loc_size_B = ds->lsize() - loc_size_A;
            PetscInt *b_indices;
            PetscMalloc(sizeof(PetscInt)*loc_size_B, &b_indices);
            for(uint i=0, j=0; i < ds->lsize(); i++)
                if (! a_used[i]) b_indices[j++] = i + ds->begin();
            ISCreateGeneral(PETSC_COMM_WORLD, loc_size_B, b_indices, PETSC_COPY_VALUES, &IsB);
            PetscFree(b_indices);
        } else {
            ISGetLocalSize(IsB, &loc_size_B);
        }
}


SchurComplement::SchurComplement(SchurComplement &other)
: LinSys_PETSC(other),
  loc_size_A(other.loc_size_A), loc_size_B(other.loc_size_B), state(other.state),
  Compl(other.Compl), ds_(other.ds_)
{
	MatCopy(other.A, A, DIFFERENT_NONZERO_PATTERN);
	MatCopy(other.IA, IA, DIFFERENT_NONZERO_PATTERN);
	MatCopy(other.IAB, IAB, DIFFERENT_NONZERO_PATTERN);
	ISCopy(other.IsA, IsA);
	ISCopy(other.IsB, IsB);
	VecCopy(other.RHS1, RHS1);
	VecCopy(other.RHS2, RHS2);
	VecCopy(other.Sol1, Sol1);
	VecCopy(other.Sol2, Sol2);

	B   = NULL;
	Bt  = NULL;
	xA  = NULL;
	C   = NULL;
}



void SchurComplement::set_from_input(const Input::Record in_rec)
{
    LinSys_PETSC::set_from_input( in_rec );

    ASSERT_PTR(Compl).error();
    Compl->set_from_input( in_rec );
}

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
    START_TIMER("form schur complement");


    PetscErrorCode ierr = 0;
    MatReuse mat_reuse;        // reuse structures after first computation of schur
    MatStructure mat_subset_pattern;
    PetscScalar *sol_array;

    mat_reuse=MAT_REUSE_MATRIX;
    mat_subset_pattern=SUBSET_NONZERO_PATTERN;
    if (state==created) {
    	mat_reuse=MAT_INITIAL_MATRIX; // indicate first construction
    	mat_subset_pattern=DIFFERENT_NONZERO_PATTERN;

        // create complement system
        // TODO: introduce LS as true object, clarify its internal states
        // create RHS sub vecs RHS1, RHS2
    	VecCreateMPI(PETSC_COMM_WORLD, loc_size_A, PETSC_DETERMINE, &(RHS1));
    	VecCreateMPI(PETSC_COMM_WORLD, loc_size_B, PETSC_DETERMINE, &(RHS2));
        VecScatterCreate(rhs_, IsA, RHS1, PETSC_NULL, &rhs1sc);
        VecScatterCreate(rhs_, IsB, RHS2, PETSC_NULL, &rhs2sc);

        // create Solution sub vecs Sol1, Sol2, Compl->solution
    	VecCreateMPI(PETSC_COMM_WORLD, loc_size_A, PETSC_DETERMINE, &(Sol1));
    	VecCreateMPI(PETSC_COMM_WORLD, loc_size_B, PETSC_DETERMINE, &(Sol2));
        VecScatterCreate(solution_, IsA, Sol1, PETSC_NULL, &sol1sc);
        VecScatterCreate(solution_, IsB, Sol2, PETSC_NULL, &sol2sc);

        VecGetArray( Sol2, &sol_array );
        Compl->set_solution( sol_array );
        VecRestoreArray( Sol2, &sol_array );


    }
    DebugOut() << print_var(mat_reuse) << print_var(matrix_changed_) << print_var(state);

    // compose Schur complement
    // Petsc need some fill estimate for results of multiplication in form nnz(A*B)/(nnz(A)+nnz(B))
    // for the first Schur compl: IA*B is bounded by ( d*(d+1) )/( d*d+2*d ) <= 5/6 for d<=4
    //                            B'*IA*B bounded by ( (d+1)*(d+1) )/ ( d*(d+1) + d ) ~ 1
    // for the second Schur :      IA*B have fill ratio ~ 1.
    //                            B'*IA*B  ...         ( N/2 *(2*N-1) )/( 2 + 2*N ) <= 1.4
    // nevertheless Petsc does not allows fill ratio below 1. so we use 1.1 for the first
    // and 1.5 for the second multiplication

    // TODO:
    // In order to let PETSC allocate structure of the complement we can not perform MatGetSubMatrix
    // and MatAXPY on the same complement matrix, since one operation change the matrix structure for the other.
    //
    // Probably no way to make this optimal using high level methods. We should have our own
    // format for schur complement matrix, store local systems and perform elimination localy.
    // Or even better assembly the complement directly. (not compatible with raw P0 method)

    if (matrix_changed_) {
       	create_inversion_matrix();

       	// compute IAB=IA*B, loc_size_B removed
		ierr+=MatGetSubMatrix(matrix_, IsA, IsB, mat_reuse, &B);
		ierr+=MatMatMult(IA, B, mat_reuse, 1.0 ,&(IAB)); // 6/7 - fill estimate
		// compute xA=Bt* IAB = Bt * IA * B, locSizeA removed
		ierr+=MatGetSubMatrix(matrix_, IsB, IsA, mat_reuse, &(Bt));
		ierr+=MatMatMult(Bt, IAB, mat_reuse, 1.9 ,&(xA)); // 1.1 - fill estimate (PETSC report values over 1.8)

		// get C block, loc_size_B removed
		ierr+=MatGetSubMatrix( matrix_, IsB, IsB, mat_reuse, &C);

		if (state==created) MatDuplicate(C, MAT_DO_NOT_COPY_VALUES, const_cast<Mat *>( Compl->get_matrix() ) );
		MatZeroEntries( *( Compl->get_matrix()) );

		// compute complement = (-1)cA+xA = Bt*IA*B - C
		if ( is_negative_definite() ) {
		    ierr+=MatAXPY(*( Compl->get_matrix() ), 1, C, SUBSET_NONZERO_PATTERN);
			ierr+=MatAXPY(*( Compl->get_matrix() ), -1, xA, mat_subset_pattern);
		} else {
			ierr+=MatAXPY(*( Compl->get_matrix() ), -1, C, SUBSET_NONZERO_PATTERN);
			ierr+=MatAXPY(*( Compl->get_matrix() ), 1, xA, mat_subset_pattern);
		}
		Compl->set_matrix_changed();

		OLD_ASSERT( ierr == 0, "PETSC Error during calculation of Schur complement.\n");

    }

    form_rhs();

	matrix_changed_ = false;

    state=formed;
}

void SchurComplement::form_rhs()
{
    START_TIMER("form rhs");
	if (rhs_changed_ || matrix_changed_) {
	    VecScatterBegin(rhs1sc, rhs_, RHS1, INSERT_VALUES, SCATTER_FORWARD);
	    VecScatterEnd(  rhs1sc, rhs_, RHS1, INSERT_VALUES, SCATTER_FORWARD);
	    VecScatterBegin(rhs2sc, rhs_, RHS2, INSERT_VALUES, SCATTER_FORWARD);
	    VecScatterEnd(  rhs2sc, rhs_, RHS2, INSERT_VALUES, SCATTER_FORWARD);

	    MatMultTranspose(IAB, RHS1, *( Compl->get_rhs() ));
	    VecAXPY(*( Compl->get_rhs() ), -1, RHS2);
	    if ( is_negative_definite() ) {
	    	VecScale(*( Compl->get_rhs() ), -1.0);
	    }
	    Compl->set_rhs_changed();
	    rhs_changed_ = false;
	}

    state=formed;
}



void SchurComplement::set_tolerances(double  r_tol, double a_tol, unsigned int max_it)
{
    LinSys_PETSC::set_tolerances(r_tol, a_tol, max_it);
    if (Compl !=nullptr) Compl->set_tolerances(r_tol, a_tol, max_it);
}

void SchurComplement::set_complement(LinSys_PETSC *ls)
{
	ASSERT_PTR(ls).error();
    Compl = ls;
}


Distribution *SchurComplement::make_complement_distribution()
{
	if (ds_ != NULL) delete ds_;
    ds_ = new Distribution(loc_size_B, PETSC_COMM_WORLD);
	return ds_;
}

void SchurComplement::create_inversion_matrix()
{
    START_TIMER("create inversion matrix");
    PetscInt ncols, pos_start, pos_start_IA;

    MatReuse mat_reuse=MAT_REUSE_MATRIX;
    if (state==created) mat_reuse=MAT_INITIAL_MATRIX; // indicate first construction

    MatGetSubMatrix(matrix_, IsA, IsA, mat_reuse, &A);
    MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &IA);
    //MatGetSubMatrix(matrix_, IsA, IsA, mat_reuse, &IA);

    MatGetOwnershipRange(A,&pos_start,PETSC_NULL);
    MatGetOwnershipRange(IA,&pos_start_IA,PETSC_NULL);

    std::vector<PetscInt> submat_rows;
    const PetscInt *cols;
    const PetscScalar *vals;

    std::vector<unsigned int> processed_rows(loc_size_A,0);

    unsigned int mat_block=1;   //actual processed block of matrix
    for(unsigned int loc_row=0; loc_row < processed_rows.size(); loc_row++) {
        if (processed_rows[loc_row] != 0) continue;

        PetscInt min=std::numeric_limits<int>::max(), max=-1, size_submat;
        PetscInt b_vals = 0; // count of values stored in B-block of Orig system
        submat_rows.clear();
        MatGetRow(A, loc_row + pos_start, &ncols, &cols, PETSC_NULL);
        for (PetscInt i=0; i<ncols; i++) {
            if (cols[i] < pos_start || cols[i] >= pos_start+loc_size_A) {
                b_vals++;
            } else {
                if (cols[i] < min) {
                    min=cols[i];
                }
                if (cols[i] > max) {
                    max=cols[i];
                }
            }
        }
        size_submat = max - min + 1;
        OLD_ASSERT(ncols-b_vals == size_submat, "Submatrix cannot contains empty values.\n");

        MatRestoreRow(A, loc_row + pos_start, &ncols, &cols, PETSC_NULL);
        arma::mat submat2(size_submat, size_submat);
        submat2.zeros();
        for (PetscInt i=0; i<size_submat; i++) {
            processed_rows[ loc_row + i ] = mat_block;
            submat_rows.push_back( i + loc_row + pos_start_IA );
            MatGetRow(A, i + loc_row + pos_start, &ncols, &cols, &vals);
            for (PetscInt j=0; j<ncols; j++) {
                if (cols[j] >= pos_start && cols[j] < pos_start+loc_size_A) {
                    submat2( i, cols[j] - loc_row - pos_start ) = vals[j];
                }
            }
            MatRestoreRow(A, i + loc_row + pos_start, &ncols, &cols, &vals);
		}
        // get inversion matrix
        arma::mat invmat = submat2.i();
        // stored to inversion IA matrix
        const PetscInt* rows = &submat_rows[0];
        MatSetValues(IA, submat_rows.size(), rows, submat_rows.size(), rows, invmat.memptr(), INSERT_VALUES);

        mat_block++;
    }

    MatAssemblyBegin(IA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(IA, MAT_FINAL_ASSEMBLY);


    /*
    MatCreateMPIAIJSumSeqAIJ(PETSC_COMM_WORLD, locIA, PETSC_DECIDE, PETSC_DECIDE, reuse, &IA);
    */
}


double SchurComplement::get_solution_precision()
{
	if (Compl != NULL) {
		return Compl->get_solution_precision();
	}
	return std::numeric_limits<double>::infinity();
}


LinSys::SolveInfo SchurComplement::solve() {
    START_TIMER("SchurComplement::solve");
    this->form_schur();

    VecScatterBegin(sol1sc, solution_, Sol1, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(  sol1sc, solution_, Sol1, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBegin(sol2sc, solution_, Sol2, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(  sol2sc, solution_, Sol2, INSERT_VALUES, SCATTER_FORWARD);
    
    //output schur complement in matlab file
//     string output_file = FilePath("schur.m", FilePath::output_file);
//     PetscViewer    viewer;
//     PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_file.c_str(), &viewer);
//     PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//     MatView( *const_cast<Mat*>(Compl->get_matrix()), viewer);
//     VecView( *const_cast<Vec*>(Compl->get_rhs()), viewer);
    
    LinSys::SolveInfo si = Compl->solve();
//     VecView(Compl->get_solution(), viewer);

    // TODO: Resolve step is not necessary inside of nonlinear solver. Can optimize.
    this->resolve();
    return si;
}


/**
 * COMPUTE ELIMINATED PART OF THE ORIG. SYS. & RESTORE RHS and SOLUTION VECTORS
 *  x_1 = IA * RHS_1 - IAB * x_2
 */

void SchurComplement::resolve()
{
    this->form_schur();

    START_TIMER("SchurComplemet::resolve without form schur");

    chkerr(MatMult(IAB,Compl->get_solution(),Sol1));
    chkerr(VecScale(Sol1,-1));
    chkerr(MatMultAdd(IA,RHS1,Sol1,Sol1));
    VecScatterBegin(sol1sc, Sol1, solution_, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(  sol1sc, Sol1, solution_, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterBegin(sol2sc, Sol2, solution_, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(  sol2sc, Sol2, solution_, INSERT_VALUES, SCATTER_REVERSE);
}


double SchurComplement::compute_residual()
{
    //DebugOut() << print_var(LinSys_PETSC::compute_residual());
    resolve();
    return LinSys_PETSC::compute_residual();
}



/**
 * SCHUR COMPLEMENT destructor
 */
SchurComplement :: ~SchurComplement() {

    if ( A  != NULL )             chkerr(MatDestroy(&A));
    if ( B  != NULL )             chkerr(MatDestroy(&B));
    if ( Bt != NULL )             chkerr(MatDestroy(&Bt));
    if ( C != NULL )              chkerr(MatDestroy(&C));
    if ( xA != NULL )             chkerr(MatDestroy(&xA));
    if ( IA != NULL )             chkerr(MatDestroy(&IA));
    if ( IAB != NULL )            chkerr(MatDestroy(&IAB));
    if ( IsA != NULL )            chkerr(ISDestroy(&IsA));
    if ( IsB != NULL )            chkerr(ISDestroy(&IsB));
    if ( RHS1 != NULL )           chkerr(VecDestroy(&RHS1));
    if ( RHS2 != NULL )           chkerr(VecDestroy(&RHS2));
    if ( Sol1 != NULL )           chkerr(VecDestroy(&Sol1));
    if ( Sol2 != NULL )           chkerr(VecDestroy(&Sol2));
    if ( rhs1sc != NULL )         chkerr(VecScatterDestroy(&rhs1sc));
    if ( rhs2sc != NULL )         chkerr(VecScatterDestroy(&rhs2sc));
    if ( sol1sc != NULL )         chkerr(VecScatterDestroy(&sol1sc));
    if ( sol2sc != NULL )         chkerr(VecScatterDestroy(&sol2sc));
    if ( IA != NULL )             chkerr(MatDestroy(&IA));

    if (Compl != NULL)            delete Compl;
    if (ds_ != NULL)              delete ds_;

}
