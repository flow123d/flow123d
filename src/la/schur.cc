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

#include "la/distribution.hh"
#include "la/local_to_global_map.hh"
#include "system/system.hh"
#include "la/linsys.hh"
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

SchurComplement::SchurComplement(IS ia, Distribution *ds)
: LinSys_PETSC(ds), IsA(ia), state(created)
{
        // check index set
        ASSERT(IsA != NULL, "Index set IsA is not defined.\n" );

        // initialize variables
        Compl   = NULL;
        IA      = NULL;
        B       = NULL;
        Bt      = NULL;
        xA      = NULL;
        IAB     = NULL;
        IsB     = NULL;
        RHS1    = NULL;
        RHS2    = NULL;
        Sol1    = NULL;
        Sol2    = NULL;

        // create A block index set
        ISGetLocalSize(IsA, &loc_size_A);

        // create B block index set
        loc_size_B = rows_ds_->lsize() - loc_size_A;
        ISCreateStride(PETSC_COMM_WORLD,loc_size_B,rows_ds_->begin()+loc_size_A,1,&IsB);
}


SchurComplement::SchurComplement(SchurComplement &other)
: LinSys_PETSC(other),
  loc_size_A(other.loc_size_A), loc_size_B(other.loc_size_B), state(other.state),
  Compl(other.Compl), ds_(other.ds_)
{
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
    PetscErrorCode ierr = 0;
    MatReuse mat_reuse;        // reuse structures after first computation of schur
    PetscScalar *rhs_array, *sol_array;

    mat_reuse=MAT_REUSE_MATRIX;
    if (state==created) {
    	mat_reuse=MAT_INITIAL_MATRIX; // indicate first construction

        // create complement system
        // TODO: introduce LS as true object, clarify its internal states
        // create RHS sub vecs RHS1, RHS2
        VecGetArray(rhs_, &rhs_array);
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,rhs_array,&(RHS1));

        // create Solution sub vecs Sol1, Compl->solution
        VecGetArray(solution_, &sol_array);
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_A,PETSC_DETERMINE,sol_array,&(Sol1));

        VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_B,PETSC_DETERMINE,rhs_array+loc_size_A,&(RHS2));
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1,loc_size_B,PETSC_DETERMINE,sol_array+loc_size_A,&(Sol2));

        VecRestoreArray(rhs_, &rhs_array);
        VecRestoreArray(solution_, &sol_array);

        VecGetArray( Sol2, &sol_array );
        Compl->set_solution( sol_array );
        Compl->set_from_input( in_rec_ );
        VecRestoreArray( Sol2, &sol_array );

    }

    // compose Schur complement
    // Petsc need some fill estimate for results of multiplication in form nnz(A*B)/(nnz(A)+nnz(B))
    // for the first Schur compl: IA*B is bounded by ( d*(d+1) )/( d*d+2*d ) <= 5/6 for d<=4
    //                            B'*IA*B bounded by ( (d+1)*(d+1) )/ ( d*(d+1) + d ) ~ 1
    // for the second Schur :      IA*B have fill ratio ~ 1.
    //                            B'*IA*B  ...         ( N/2 *(2*N-1) )/( 2 + 2*N ) <= 1.4
    // nevertheless Petsc does not allows fill ratio below 1. so we use 1.1 for the first
    // and 1.5 for the second multiplication

    if (matrix_changed_) {
       	create_inversion_matrix();

       	// compute IAB=IA*B, loc_size_B removed
		ierr+=MatGetSubMatrix(matrix_, IsA, IsB, mat_reuse, &B);
		ierr+=MatMatMult(IA, B, mat_reuse, 1.0 ,&(IAB)); // 6/7 - fill estimate
		// compute xA=Bt* IAB = Bt * IA * B, locSizeA removed
		ierr+=MatGetSubMatrix(matrix_, IsB, IsA, mat_reuse, &(Bt));
		ierr+=MatMatMult(Bt, IAB, mat_reuse, 1.9 ,&(xA)); // 1.1 - fill estimate (PETSC report values over 1.8)

		// get C block, loc_size_B removed
		ierr+=MatGetSubMatrix( matrix_, IsB, IsB, mat_reuse, const_cast<Mat *>( Compl->get_matrix() ) );
		// compute complement = (-1)cA+xA = Bt*IA*B - C
		if ( is_negative_definite() ) {
			ierr+=MatAXPY(*( Compl->get_matrix() ), -1, xA, SUBSET_NONZERO_PATTERN);
		} else {
			ierr+=MatScale(*( Compl->get_matrix() ),-1.0);
			ierr+=MatAXPY(*( Compl->get_matrix() ), 1, xA, SUBSET_NONZERO_PATTERN);
		}
		Compl->set_matrix_changed();

		ASSERT( ierr == 0, "PETSC Error during calculation of Schur complement.\n");

    }

    form_rhs();

	matrix_changed_ = false;

    state=formed;
}

void SchurComplement::form_rhs()
{
	if (rhs_changed_ || matrix_changed_) {
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


/**
 * COMPUTE ELIMINATED PART OF THE ORIG. SYS. & RESTORE RHS and SOLUTION VECTORS
 *  x_1 = IA * RHS_1 - IAB * x_2
 */

void SchurComplement::resolve()
{
    MatMult(IAB,Compl->get_solution(),Sol1);

    VecScale(Sol1,-1);

    MatMultAdd(IA,RHS1,Sol1,Sol1);

}

void SchurComplement::set_complement(LinSys_PETSC *ls)
{
	ASSERT(ls != nullptr, "NULL complement ls.\n");
    Compl = ls;
}


Distribution *SchurComplement::make_complement_distribution()
{
    ds_ = new Distribution(loc_size_B, PETSC_COMM_WORLD);
	return ds_;
}

void SchurComplement::create_inversion_matrix()
{
    PetscInt ncols, pos_start, pos_start_IA;

    MatReuse mat_reuse=MAT_REUSE_MATRIX;
    if (state==created) mat_reuse=MAT_INITIAL_MATRIX; // indicate first construction

    MatGetSubMatrix(matrix_, IsA, IsA, mat_reuse, &IA);
    MatGetOwnershipRange(matrix_,&pos_start,PETSC_NULL);
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
        MatGetRow(matrix_, loc_row + pos_start, &ncols, &cols, PETSC_NULL);
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
        ASSERT(ncols-b_vals == size_submat, "Submatrix cannot contains empty values.\n");

        MatRestoreRow(matrix_, loc_row + pos_start, &ncols, &cols, PETSC_NULL);
        arma::mat submat2(size_submat, size_submat);
        submat2.zeros();
        for (PetscInt i=0; i<size_submat; i++) {
            processed_rows[ loc_row + i ] = mat_block;
            submat_rows.push_back( i + loc_row + pos_start_IA );
            MatGetRow(matrix_, i + loc_row + pos_start, &ncols, &cols, &vals);
            for (PetscInt j=0; j<ncols; j++) {
                if (cols[j] >= pos_start && cols[j] < pos_start+loc_size_A) {
                    submat2( i, cols[j] - loc_row - pos_start ) = vals[j];
                }
            }
            MatRestoreRow(matrix_, i + loc_row + pos_start, &ncols, &cols, &vals);
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
}


double SchurComplement::get_solution_precision()
{
	if (Compl != NULL) {
		return Compl->get_solution_precision();
	}
	return std::numeric_limits<double>::infinity();
}


int SchurComplement::solve() {
	this->form_schur();

	int converged_reason = Compl->solve();
	this->resolve();

	return converged_reason;
}


/**
 * SCHUR COMPLEMENT destructor
 */
SchurComplement :: ~SchurComplement() {

    if ( B  != NULL )             MatDestroy(&B);
    if ( Bt != NULL )             MatDestroy(&Bt);
    if ( xA != NULL )             MatDestroy(&xA);
    if ( IA != NULL )             MatDestroy(&IA);
    if ( IAB != NULL )            MatDestroy(&IAB);
    if ( IsA != NULL )            ISDestroy(&IsA);
    if ( IsB != NULL )            ISDestroy(&IsB);
    if ( RHS1 != NULL )           VecDestroy(&RHS1);
    if ( RHS2 != NULL )           VecDestroy(&RHS2);
    if ( Sol1 != NULL )           VecDestroy(&Sol1);
    if ( Sol2 != NULL )           VecDestroy(&Sol2);
    if ( IA != NULL )             MatDestroy(&IA);

    if (Compl != NULL)            delete Compl;

}
