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
 * @file    linsys_PETSC.cc
 * @brief   Solver based on the original PETSc solver using MPIAIJ matrix and succesive Schur complement construction
 * @author  Jakub Sistek
 */

// derived from base linsys
#include "la/linsys_PERMON.hh"
#include "petscvec.h"
#include "petscksp.h"
#include "petscmat.h"
#include "system/sys_profiler.hh"
#include "system/system.hh"
#include "fem/dofhandler.hh"



namespace it = Input::Type;

const it::Record & LinSys_PERMON::get_input_type() {
	return it::Record("Permon", "PERMON solver settings.\n It provides interface to various PERMON solvers. The convergence criteria is:\n"
	        "```\n"
	        "norm( res_i )  < max( norm( res_0 ) * r_tol, a_tol )\n"
	        "```\n"
	        "where ```res_i``` is the residuum vector after i-th iteration of the solver and ```res_0``` is the estimate of the norm of the initial residual. "
	        "If the initial guess of the solution is provided (usually only for transient equations) the residual of this estimate is used, "
	        "otherwise the norm of preconditioned RHS is used. "
	        "The default norm is (($L_2$)) norm of preconditioned residual: (($ P^{-1}(Ax-b)$)), usage of other norm may be prescribed using the 'option' key. "
	        "See also PETSc documentation for KSPSetNormType.")
		.derive_from(LinSys::get_input_type())
		.declare_key("r_tol", it::Double(0.0, 1.0), it::Default::read_time("Default value is set by the nonlinear solver or the equation. "
                        "If not, we use the value 1.0e-7."),
					"Residual tolerance relative to the initial error.")
		.declare_key("a_tol", it::Double(0.0), it::Default::read_time("Default value is set by the nonlinear solver or the equation. "
                        "If not, we use the value 1.0e-11."),
		            "Absolute residual tolerance.")
        .declare_key("d_tol", it::Double(0.0), it::Default::read_time("Default value is set by the nonlinear solver or the equation. "
                        "If not, we use the value 10000."),
		            "Tolerance for divergence.")
        .declare_key("max_it", it::Integer(0), it::Default::read_time("Default value is set by the nonlinear solver or the equation. "
                        "If not, we use the value 1000."),
                    "Maximum number of outer iterations of the linear solver.")
    .declare_key("warm_start", it::Bool(), it::Default("true"), "Warm start QPS solver with the privous solution.")
		.declare_key("options", it::String(), it::Default("\"\""),  "This options is passed to PETSC to create a particular KSP (Krylov space method).\n"
                                                                    "If the string is left empty (by default), the internal default options is used.")
		.close();
}


const int LinSys_PERMON::registrar = LinSys_PERMON::get_input_type().size();


LinSys_PERMON::LinSys_PERMON(const DOFHandlerMultiDim &dh, const std::string &params)
        : LinSys( dh.distr().get() ),
          params_(params),
          init_guess_nonzero(false),
          matrix_(NULL),
          use_feti_(false)
{
    PetscErrorCode ierr;

    // create PETSC vector for rhs
    ierr = VecCreateMPI( comm_, rows_ds_->lsize(), PETSC_DECIDE, &rhs_ ); CHKERRV( ierr );
    ierr = VecZeroEntries( rhs_ ); CHKERRV( ierr );
    VecDuplicate(rhs_, &residual_);

    // create l2g map
    ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, dh.get_local_to_global_map().size(), dh.get_local_to_global_map().data(), PETSC_USE_POINTER, &l2g_); CHKERRV(ierr);

    solution_precision_ = std::numeric_limits<double>::infinity();
    matrix_changed_ = true;
    rhs_changed_ = true;

    matrix_ineq_ = NULL;
    ineq_ = NULL;
    warm_solution_ = NULL;
    maxeig_ = PETSC_DECIDE;
}

LinSys_PERMON::LinSys_PERMON( LinSys_PERMON &other )
	: LinSys(other), params_(other.params_), l2g_(other.l2g_), solution_precision_(other.solution_precision_)
{
    MatCopy(other.matrix_, matrix_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.rhs_, rhs_);
	VecCopy(other.on_vec_, on_vec_);
	VecCopy(other.off_vec_, off_vec_);

	MatCopy(other.matrix_ineq_, matrix_ineq_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.ineq_, ineq_);
	VecCopy(other.warm_solution_, warm_solution_);
    warm_start_ = other.warm_start_;
    maxeig_ = other.maxeig_;
    use_feti_ = other.use_feti_;
}


void LinSys_PERMON::set_inequality(Mat matrix_ineq, Vec ineq)
{
  // TODO ref count?
  matrix_ineq_ = matrix_ineq;
  ineq_ = ineq;
}


void LinSys_PERMON::set_tolerances(double  r_tol, double a_tol, double d_tol, unsigned int max_it)
{
    if (! in_rec_.is_empty()) {
        // input record is set
        r_tol_ = in_rec_.val<double>("r_tol", r_tol);
        a_tol_ = in_rec_.val<double>("a_tol", a_tol);
        d_tol_ = in_rec_.val<double>("d_tol", d_tol);
        max_it_ = in_rec_.val<unsigned int>("max_it", max_it);
    } else {
        r_tol_ = r_tol;
        a_tol_ = a_tol;
        d_tol_ = d_tol;
        max_it_ = max_it;

    }
}


void LinSys_PERMON::start_allocation( )
{
    PetscErrorCode ierr;

    ierr = VecCreateMPI( comm_, rows_ds_->lsize(), PETSC_DECIDE, &(on_vec_) ); CHKERRV( ierr ); 
    ierr = VecSetLocalToGlobalMapping( on_vec_, l2g_ ); CHKERRV( ierr );
    ierr = VecDuplicate( on_vec_, &(off_vec_) ); CHKERRV( ierr ); 
    status_ = ALLOCATE;
}

void LinSys_PERMON::start_add_assembly()
{
    switch ( status_ ) {
        case ALLOCATE:
            this->preallocate_matrix( );
            break;
        case INSERT:
            this->finish_assembly( MAT_FLUSH_ASSEMBLY );
            break;
        case ADD:
        case DONE:
            break;
        default:
        	ASSERT_PERMANENT(false).error("Can not set values. Matrix is not preallocated.\n");
    }
    status_ = ADD;
}

void LinSys_PERMON::start_insert_assembly()
{
    switch ( status_ ) {
        case ALLOCATE:
            this->preallocate_matrix();
            break;
        case ADD:
            this->finish_assembly( MAT_FLUSH_ASSEMBLY );
            break;
        case INSERT:
        case DONE:
            break;
        default:
        	ASSERT_PERMANENT(false).error("Can not set values. Matrix is not preallocated.\n");
    }
    status_ = INSERT;
}


void LinSys_PERMON::mat_set_values_local( int nrow, int *rows, int ncol, int *cols, double *vals )
{
    // here vals would need to be converted from double to PetscScalar if it was ever something else than double :-)
    switch (status_) {
        case INSERT:
        case ADD:
            chkerr(MatSetValuesLocal(matrix_,nrow,rows,ncol,cols,vals,(InsertMode)status_));
            break;
        case ALLOCATE:
            this->preallocate_values_local(nrow,rows,ncol,cols); 
            break;
        default: DebugOut() << "LS SetValues with non allowed insert mode.\n";
    }

    matrix_changed_ = true;
}


void LinSys_PERMON::rhs_set_values( int nrow, int *rows, double *vals )
{
    PetscErrorCode ierr;

    switch (status_) {
        case INSERT:
        case ADD:
            ierr = VecSetValues(rhs_,nrow,rows,vals,(InsertMode)status_); CHKERRV( ierr ); 
            break;
        case ALLOCATE: 
            break;
        default: ASSERT_PERMANENT(false).error("LinSys's status disallow set values.\n");
    }

    rhs_changed_ = true;
}


void LinSys_PERMON::preallocate_values_local(int nrow,int *rows,int ncol,int *)
{
    for (int i=0; i<nrow; i++)
        VecSetValueLocal(on_vec_,rows[i],(double)ncol,ADD_VALUES);
}


void LinSys_PERMON::finish_assembly( )
{
    MatAssemblyType assemblyType = MAT_FINAL_ASSEMBLY;
    this->finish_assembly( assemblyType );
}


void LinSys_PERMON::finish_assembly( MatAssemblyType assembly_type )
{
    PetscErrorCode ierr;

    if (status_ == ALLOCATE) {
    	WarningOut() << "Finalizing linear system without setting values.\n";
        this->preallocate_matrix();
    }
    ierr = MatAssemblyBegin(matrix_, assembly_type); CHKERRV( ierr ); 
    ierr = VecAssemblyBegin(rhs_); CHKERRV( ierr ); 
    ierr = MatAssemblyEnd(matrix_, assembly_type); CHKERRV( ierr ); 
    ierr = VecAssemblyEnd(rhs_); CHKERRV( ierr ); 

    if (assembly_type == MAT_FINAL_ASSEMBLY) status_ = DONE;

    matrix_changed_ = true;
    rhs_changed_ = true;
}


void LinSys_PERMON::preallocate_matrix()
{
	ASSERT_EQ(status_, ALLOCATE).error("Linear system has to be in ALLOCATE status.");

    PetscErrorCode ierr;
    PetscInt *on_nz, *off_nz;
    PetscScalar *on_array, *off_array;

    // assembly and get values from counting vectors, destroy them
    VecAssemblyBegin(on_vec_);
    VecAssemblyBegin(off_vec_);

    unsigned int lsize;
    ISLocalToGlobalMappingGetSize(l2g_, (int*)(&lsize));

    on_nz  = new PetscInt[ lsize ];
    off_nz = new PetscInt[ lsize ];

    VecAssemblyEnd(on_vec_);
    VecAssemblyEnd(off_vec_);

    VecGetArray( on_vec_,  &on_array );
    VecGetArray( off_vec_, &off_array );

    for ( unsigned int i=0; i<lsize; i++ ) {
        on_nz[i]  = std::min( lsize, static_cast<uint>( on_array[i]+0.1  ) );  // small fraction to ensure correct rounding
        off_nz[i] = std::min( rows_ds_->size() - lsize, static_cast<uint>( off_array[i]+0.1 ) );
    }

    VecRestoreArray(on_vec_,&on_array);
    VecRestoreArray(off_vec_,&off_array);
    VecDestroy(&on_vec_);
    VecDestroy(&off_vec_);

    // create PETSC matrix with preallocation
    if (matrix_ != NULL)
    {
    	chkerr(MatDestroy(&matrix_));
    }
    ierr = MatCreateIS(PETSC_COMM_WORLD, 1, rows_ds_->lsize(), rows_ds_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                                  l2g_, l2g_, &matrix_); CHKERRV( ierr );
    ierr = MatISSetPreallocation(matrix_, 0, on_nz, 0, off_nz);

    if (symmetric_) MatSetOption(matrix_, MAT_SYMMETRIC, PETSC_TRUE);
    MatSetOption(matrix_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    // This option is used in order to assembly larger local matrices with own non-zero structure.
    // Zero entries are ignored so we must prevent adding exact zeroes.
    // Add LocalSystem::almost_zero for entries that should not be eliminated.
    MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);



    delete[] on_nz;
    delete[] off_nz;
}


LinSys::SolveInfo LinSys_PERMON::solve()
{
    const char *petsc_dflt_opt;
    int nits;
    
    // -mat_no_inode ... inodes are usefull only for
    //  vector problems e.g. MH without Schur complement reduction
    
    /* Comment to PETSc options:
     * 
     * -ksp_diagonal_scale scales the matrix before solution, while -ksp_diagonal_scale_fix just fixes the scaling after solution
     * -pc_asm_type basic enforces classical Schwartz method, which seems more stable for positive definite systems.
     *                    The default 'restricted' probably violates s.p.d. structure, many tests fail.
     */
    if (rows_ds_->np() > 1) {
        // parallel setting
       if (this->is_positive_definite())
           petsc_dflt_opt="-ksp_type cg -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_type asm -pc_asm_type basic -pc_asm_overlap 4 -sub_pc_type icc -sub_pc_factor_levels 3  -sub_pc_factor_fill 6.0";
           //petsc_dflt_opt="-ksp_type bcgs -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3  -sub_pc_factor_fill 6.0";
       else
           petsc_dflt_opt="-ksp_type bcgs -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3 -sub_pc_factor_fill 6.0";
    
    } 
    else {
        // serial setting
       if (this->is_positive_definite())
           petsc_dflt_opt="-ksp_type cg -pc_type icc  -pc_factor_levels 3 -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    	   //petsc_dflt_opt="-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
       else
           petsc_dflt_opt="-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    }

    if (params_ == "") params_ = petsc_dflt_opt;
    LogOut().fmt("inserting petsc options: {}\n",params_.c_str());
    
    // now takes an optional PetscOptions object as the first argument
    // value NULL will preserve previous behaviour previous behavior.
    PetscOptionsInsertString(NULL, params_.c_str()); // overwrites previous options values
    
    MatSetOption( matrix_, MAT_USE_INODES, PETSC_FALSE );

    if (use_feti_ && (rows_ds_->np() == 1)) {
        WarningOut() << "FETI can be only used on multiple processes - switching to standard QP solver.";
        use_feti_ = false;
    }

    
    chkerr(QPCreate(comm_, &system));
    chkerr(QPSetOptionsPrefix(system,"permon_")); // avoid clash on PC objects from hydro PETSc solver
    if (use_feti_) {

        // convert to MatIS format
        Mat Afixed;
        chkerr(MatConvert(matrix_, MATIS, MAT_INITIAL_MATRIX, &Afixed));

        // create new matrix without zero rows
        Mat Aloc, Aloc_feti;
        IS isnz, isnz_tmp;
        ISLocalToGlobalMapping l2g, mapping;
        const PetscInt *l2g_arr, *nz_arr;
        PetscInt *l2g_feti_arr;
        PetscInt nmap, nrows, i, j, k;
        
        chkerr(MatISGetLocalMat(Afixed, &Aloc));
        chkerr(MatFindNonzeroRows(Aloc, &isnz_tmp));
        if (isnz_tmp != NULL) {
            const PetscInt *nz_idx;
            ISGetLocalSize(isnz_tmp, &nrows);
            ISGetIndices(isnz_tmp, &nz_idx);
            ISCreateGeneral(PETSC_COMM_WORLD, nrows, nz_idx, PETSC_COPY_VALUES, &isnz);
            ISRestoreIndices(isnz_tmp, &nz_idx);
        } else {
            MatGetLocalSize(Aloc, &nrows, NULL);
            PetscInt *nz_idx;
            PetscMalloc1(nrows, &nz_idx);
            for (int i=0; i<nrows; i++)
                nz_idx[i] = i;
            ISCreateGeneral(PETSC_COMM_WORLD, nrows, nz_idx, PETSC_OWN_POINTER, &isnz);
        }
        ISDestroy(&isnz_tmp);
        chkerr(MatCreateSubMatrix(Aloc, isnz, isnz, MAT_INITIAL_MATRIX, &Aloc_feti));
        chkerr(MatISRestoreLocalMat(Afixed, &Aloc));
        chkerr(ISGetIndices(isnz, &nz_arr));
        chkerr(ISGetLocalSize(isnz, &nrows));
        chkerr(MatISGetLocalToGlobalMapping(Afixed, &l2g, NULL));
        chkerr(ISLocalToGlobalMappingGetIndices(l2g, &l2g_arr));
        chkerr(ISLocalToGlobalMappingGetSize(l2g, &nmap));
        chkerr(PetscMalloc1(nrows, &l2g_feti_arr));
        k = 0;
        for (i=0; i<=nmap; i++) {
            for (j=0; j<=nrows; j++) {
                if (i == nz_arr[j]) {
                    l2g_feti_arr[k] = l2g_arr[i];
                    k++;
                    break;
                }
            }
        }
        chkerr(ISRestoreIndices(isnz, &nz_arr));
        chkerr(ISDestroy(&isnz));
        chkerr(ISLocalToGlobalMappingCreate(comm_, 1, nrows, l2g_feti_arr, PETSC_OWN_POINTER, &mapping));
        chkerr(ISLocalToGlobalMappingRestoreIndices(l2g, &l2g_arr));
        chkerr(MatDestroy(&Afixed));
        chkerr(MatCreateIS(PETSC_COMM_WORLD, 1, rows_ds_->lsize(), rows_ds_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                            mapping, mapping, &Afixed));
        chkerr(MatISSetLocalMat(Afixed, Aloc_feti));
        chkerr(MatAssemblyBegin(Afixed, MAT_FINAL_ASSEMBLY));
        chkerr(MatAssemblyEnd(Afixed, MAT_FINAL_ASSEMBLY));
        chkerr(MatDestroy(&Aloc_feti));
        chkerr(ISLocalToGlobalMappingDestroy(&mapping));

        // set QP matrix
        chkerr(QPSetOperator(system, Afixed));
        chkerr(MatDestroy(&Afixed));
    } else {
        Mat matrix_aij;
        chkerr(MatConvert(matrix_, MATAIJ, MAT_INITIAL_MATRIX, &matrix_aij));
        chkerr(QPSetOperator(system, matrix_aij));
    }
    chkerr(QPSetRhs(system, rhs_));
    chkerr(QPSetInitialVector(system, solution_));
    if (ineq_) {
      chkerr(QPSetIneq(system, matrix_ineq_, ineq_));
    }
    if (use_feti_) { // FETI
      chkerr(QPTMatISToBlockDiag(system));
      chkerr(QPGetChild(system, &system));
      chkerr(QPFetiSetUp(system));
      chkerr(PetscOptionsInsertString(NULL, "-feti"));
    } else if (ineq_) { // dualization without FETI
      chkerr(QPTDualize(system, MAT_INV_MONOLITHIC, MAT_REG_NONE));
    }
      
    // Set/Unset additional transformations, e.g -project 0 for projector avoiding FETI
    if (use_feti_) {
        chkerr(QPTFromOptions(system));
        chkerr(QPGetParent(system, &system));
    }

    // Set runtime options, e.g -qp_chain_view_kkt
    chkerr(QPSetFromOptions(system));
    
    chkerr(QPSCreate(comm_, &solver));
    chkerr(QPSSetQP(solver, system));

    // TODO take care of tolerances - shall we support both input file and command line petsc setting
    chkerr(QPSSetTolerances(solver, r_tol_, a_tol_, d_tol_, PETSC_DEFAULT));
    chkerr(QPSSetTolerances(solver, r_tol_, a_tol_, d_tol_,  max_it_));
    chkerr(QPSSetFromOptions(solver));
    chkerr(QPSSetUp(solver)); // solvers may do additional QP tranformations

    // warm start the solver
    if (warm_solution_ != NULL) {
      QP qp_last;
      Vec sol;
      chkerr(QPChainGetLast(system,&qp_last));
      chkerr(QPGetSolutionVector(qp_last,&sol));
      chkerr(VecCopy(warm_solution_,sol));
    }

    KSPConvergedReason reason;
    {
		START_TIMER("PERMON linear solver");
		START_TIMER("PERMON linear iteration");
		chkerr(QPSSolve(solver));
		QPSGetConvergedReason(solver,&reason);
		QPSGetIterationNumber(solver,&nits);
		ADD_CALLS(nits);
    }
    // substitute by PETSc call for residual
    //VecNorm(rhs_, NORM_2, &residual_norm_);
    
    LogOut().fmt("convergence reason {}, number of iterations is {}\n", reason, nits);

    // get residual norm
    compute_residual();

    // TODO: I do not understand this 
    //Profiler::instance()->set_timer_subframes("SOLVING MH SYSTEM", nits);

    // set inner solver solution for warm start
    if (warm_start_) {
      QP qp_last;
      Vec sol;
      PetscBool same;
      chkerr(QPChainGetLast(system,&qp_last));
      chkerr(QPGetSolutionVector(qp_last,&sol));
      chkerr(VecDestroy(&warm_solution_));
      chkerr(VecDuplicate(sol,&warm_solution_));
      chkerr(VecCopy(sol,warm_solution_));
      chkerr(PetscObjectTypeCompare((PetscObject)solver,QPSMPGP,&same));
      if (same) {
        char stri[128];
        chkerr(QPSMPGPGetOperatorMaxEigenvalue(solver,&maxeig_));
        chkerr(PetscSNPrintf(stri, sizeof(stri), "-qps_mpgp_maxeig %.17g",maxeig_));
        // NOTE: this could be done through API, but needs to have correct
        // order for QPSSetUp/ warm start set up, that depends on the solver type
        chkerr(PetscOptionsInsertString(NULL,stri));
      } // TODO: inherit maxeig for SMALXE
    }

    chkerr(QPSDestroy(&solver));
    chkerr(QPDestroy(&system));

    return LinSys::SolveInfo(static_cast<int>(reason), static_cast<int>(nits));

}

void LinSys_PERMON::view(string text )
{
    FilePath matFileName(text + "_flow123d_matrix.m",FilePath::FileType::output_file);
    FilePath rhsFileName(text + "_flow123d_rhs.m",FilePath::FileType::output_file);
    FilePath solFileName(text + "_flow123d_sol.m",FilePath::FileType::output_file);
    FilePath mat_ineqFileName(text + "_flow123d_matrix_ineq.m",FilePath::FileType::output_file);
    FilePath ineqFileName(text + "_flow123d_ineq.m",FilePath::FileType::output_file);

    PetscViewer myViewer;

    if ( matrix_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)matFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        MatView( matrix_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the matrix of LinSys is not set.\n";

    if ( rhs_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)rhsFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( rhs_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the rhs of LinSys is not set.\n";

    if ( solution_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)solFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( solution_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the solution of LinSys is not set.\n";
    if ( matrix_ineq_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)mat_ineqFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        MatView( matrix_ineq_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the inequality matrix of LinSys is not set.\n";
    if ( ineq_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)ineqFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( ineq_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the inequality vector of LinSys is not set.\n";
}

LinSys_PERMON::~LinSys_PERMON( )
{
    if (warm_solution_ != NULL) chkerr(VecDestroy(&warm_solution_));
}

void LinSys_PERMON::set_from_input(const Input::Record in_rec)
{
    LinSys::set_from_input( in_rec );

	// PETSC specific parameters
    // If parameters are specified in input file, they are used,
    // otherwise keep settings provided in constructor of LinSys_PETSC.
    std::string user_params = in_rec.val<string>("options");
	if (user_params != "") params_ = user_params;

    // PERMON specific parameters
    warm_start_ = false;//in_rec.val<bool>("warm_start");
    warm_start_ = true;//in_rec.val<bool>("warm_start");
}

double LinSys_PERMON::get_solution_precision()
{
	return solution_precision_;
}


double LinSys_PERMON::compute_residual()
{
    // TODO return ||A*x - b + B'*lambda||?
    MatMult(matrix_, solution_, residual_);
    VecAXPY(residual_,-1.0, rhs_);
    VecNorm(residual_, NORM_2, &solution_precision_);
    return solution_precision_;
}
