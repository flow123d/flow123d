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
#include "la/linsys_PETSC.hh"
#include "petscvec.h"
#include "petscksp.h"
#include "petscmat.h"
#include "system/sys_profiler.hh"
#include "system/system.hh"
#include "fem/dofhandler.hh"


//#include <boost/bind.hpp>

namespace it = Input::Type;

const it::Record & LinSys_PETSC::get_input_type() {
	return it::Record("Petsc", "PETSc solver settings.\n It provides interface to various PETSc solvers. The convergence criteria is:\n"
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
        .declare_key("max_it", it::Integer(0), it::Default::read_time("Default value is set by the nonlinear solver or the equation. "
                        "If not, we use the value 1000."),
                    "Maximum number of outer iterations of the linear solver.")
		.declare_key("options", it::String(), it::Default("\"\""),  "This options is passed to PETSC to create a particular KSP (Krylov space method).\n"
                                                                    "If the string is left empty (by default), the internal default options is used.")
		.close();
}


const int LinSys_PETSC::registrar = LinSys_PETSC::get_input_type().size();


LinSys_PETSC::LinSys_PETSC( const Distribution * rows_ds, const std::string &params)
        : LinSys( rows_ds ),
          params_(params),
          init_guess_nonzero(false),
          l2g_(nullptr),
          matrix_(0)
{
    // create PETSC vectors:
    PetscErrorCode ierr;
    // rhs
    v_rhs_= new double[ rows_ds_->lsize() + 1 ];
    ierr = VecCreateMPIWithArray( comm_, 1, rows_ds_->lsize(), PETSC_DECIDE, v_rhs_, &rhs_ ); CHKERRV( ierr );
    ierr = VecZeroEntries( rhs_ ); CHKERRV( ierr );
    VecDuplicate(rhs_, &residual_);

    matrix_ = NULL;
    solution_precision_ = std::numeric_limits<double>::infinity();
    matrix_changed_ = true;
    rhs_changed_ = true;
}

LinSys_PETSC::LinSys_PETSC( const DOFHandlerMultiDim &dh, const std::string &params)
        : LinSys( dh.distr().get() ),
          params_(params),
          init_guess_nonzero(false),
          l2g_(std::make_shared<const std::vector<LongIdx>>(dh.get_local_to_global_map())),
          matrix_(0)
{
    // create PETSC vectors:
    PetscErrorCode ierr;
    // rhs
    v_rhs_= new double[ rows_ds_->lsize() + 1 ];
    ierr = VecCreateMPIWithArray( comm_, 1, rows_ds_->lsize(), PETSC_DECIDE, v_rhs_, &rhs_ ); CHKERRV( ierr );
    ierr = VecZeroEntries( rhs_ ); CHKERRV( ierr );
    VecDuplicate(rhs_, &residual_);

    matrix_ = NULL;
    solution_precision_ = std::numeric_limits<double>::infinity();
    matrix_changed_ = true;
    rhs_changed_ = true;
}

LinSys_PETSC::LinSys_PETSC( LinSys_PETSC &other )
	: LinSys(other), params_(other.params_), l2g_(other.l2g_), v_rhs_(NULL), solution_precision_(other.solution_precision_)
{
	MatCopy(other.matrix_, matrix_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.rhs_, rhs_);
	VecCopy(other.on_vec_, on_vec_);
	VecCopy(other.off_vec_, off_vec_);
}

void LinSys_PETSC::set_tolerances(double  r_tol, double a_tol, unsigned int max_it)
{
    if (! in_rec_.is_empty()) {
        // input record is set
        r_tol_ = in_rec_.val<double>("r_tol", r_tol);
        a_tol_ = in_rec_.val<double>("a_tol", a_tol);
        max_it_ = in_rec_.val<unsigned int>("max_it", max_it);
    } else {
        r_tol_ = r_tol;
        a_tol_ = a_tol;
        max_it_ = max_it;

    }
}


void LinSys_PETSC::start_allocation( )
{
    PetscErrorCode ierr;

    if (l2g_ == nullptr)
    {
        ierr = VecCreateMPI( comm_, rows_ds_->lsize(), PETSC_DECIDE, &(on_vec_) ); CHKERRV( ierr ); 
        ierr = VecDuplicate( on_vec_, &(off_vec_) ); CHKERRV( ierr ); 
    }
    else
    {
        ISLocalToGlobalMapping l2g_is;
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, l2g_->size(), l2g_->data(), PETSC_USE_POINTER, &l2g_is);
        ierr = VecCreateMPI( comm_, rows_ds_->lsize(), PETSC_DECIDE, &(on_vec_) ); CHKERRV( ierr ); 
        ierr = VecSetLocalToGlobalMapping( on_vec_, l2g_is ); CHKERRV( ierr );
        ierr = VecDuplicate( on_vec_, &(off_vec_) ); CHKERRV( ierr ); 
    }
    status_ = ALLOCATE;
}

void LinSys_PETSC::start_add_assembly()
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

void LinSys_PETSC::start_insert_assembly()
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

void LinSys_PETSC::mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals )
{
    // here vals would need to be converted from double to PetscScalar if it was ever something else than double :-)
    switch (status_) {
        case INSERT:
        case ADD:
            chkerr(MatSetValues(matrix_,nrow,rows,ncol,cols,vals,(InsertMode)status_));
            break;
        case ALLOCATE:
            this->preallocate_values(nrow,rows,ncol,cols); 
            break;
        default: DebugOut() << "LS SetValues with non allowed insert mode.\n";
    }

    matrix_changed_ = true;
}

void LinSys_PETSC::rhs_set_values( int nrow, int *rows, double *vals )
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

void LinSys_PETSC::preallocate_values(int nrow,int *rows,int ncol,int *cols)
{
    int i,j;
    int col;
    PetscInt row;

    for (i=0; i<nrow; i++) {
        row=rows[i];
        for(j=0; j<ncol; j++) {
            col = cols[j];
            if (rows_ds_->get_proc(row) == rows_ds_->get_proc(col))
                VecSetValue(on_vec_,row,1.0,ADD_VALUES);
            else
                VecSetValue(off_vec_,row,1.0,ADD_VALUES);
        }
    }
}

void LinSys_PETSC::mat_set_values_local( int nrow, int *rows, int ncol, int *cols, double *vals )
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

void LinSys_PETSC::rhs_set_values_local( int nrow, int *rows, double *vals )
{
    PetscErrorCode ierr;

    switch (status_) {
        case INSERT:
        case ADD:
            ierr = VecSetValuesLocal(rhs_,nrow,rows,vals,(InsertMode)status_); CHKERRV( ierr ); 
            break;
        case ALLOCATE: 
            break;
        default: ASSERT_PERMANENT(false).error("LinSys's status disallow set values.\n");
    }

    rhs_changed_ = true;
}

void LinSys_PETSC::preallocate_values_local(int nrow,int *rows,int ncol,int *)
{
    for (int i=0; i<nrow; i++)
        VecSetValueLocal(on_vec_,rows[i],(double)ncol,ADD_VALUES);
}

void LinSys_PETSC::preallocate_matrix()
{
	ASSERT_EQ(status_, ALLOCATE).error("Linear system has to be in ALLOCATE status.");

    PetscErrorCode ierr;
    PetscInt *on_nz, *off_nz;
    PetscScalar *on_array, *off_array;

    // assembly and get values from counting vectors, destroy them
    VecAssemblyBegin(on_vec_);
    VecAssemblyBegin(off_vec_);

    on_nz  = new PetscInt[ rows_ds_->lsize() ];
    off_nz = new PetscInt[ rows_ds_->lsize() ];

    VecAssemblyEnd(on_vec_);
    VecAssemblyEnd(off_vec_);

    VecGetArray( on_vec_,  &on_array );
    VecGetArray( off_vec_, &off_array );

    for ( unsigned int i=0; i<rows_ds_->lsize(); i++ ) {
        on_nz[i]  = std::min( rows_ds_->lsize(), static_cast<uint>( on_array[i]+0.1  ) );  // small fraction to ensure correct rounding
        off_nz[i] = std::min( rows_ds_->size() - rows_ds_->lsize(), static_cast<uint>( off_array[i]+0.1 ) );
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
    if (l2g_ == nullptr)
    {
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, rows_ds_->lsize(), rows_ds_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                                0, on_nz, 0, off_nz, &matrix_); CHKERRV( ierr );
    }
    else
    {
        ISLocalToGlobalMapping l2g_is;
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, l2g_->size(), l2g_->data(), PETSC_USE_POINTER, &l2g_is);
        ISLocalToGlobalMappingView( l2g_is, PETSC_VIEWER_STDOUT_WORLD );
        ierr = MatCreateIS(PETSC_COMM_WORLD, 1, rows_ds_->lsize(), rows_ds_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                                  l2g_is, l2g_is, &matrix_); CHKERRV( ierr );
        ierr = MatISSetPreallocation(matrix_, 0, on_nz, 0, off_nz);
    }

    if (symmetric_) MatSetOption(matrix_, MAT_SYMMETRIC, PETSC_TRUE);
    MatSetOption(matrix_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    // This option is used in order to assembly larger local matrices with own non-zero structure.
    // Zero entries are ignored so we must prevent adding exact zeroes.
    // Add LocalSystem::almost_zero for entries that should not be eliminated.
    MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);



    delete[] on_nz;
    delete[] off_nz;
}

void LinSys_PETSC::finish_assembly( )
{
    MatAssemblyType assemblyType = MAT_FINAL_ASSEMBLY;
    this->finish_assembly( assemblyType );
}

void LinSys_PETSC::finish_assembly( MatAssemblyType assembly_type )
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

    //PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INDEX);
    //MatView(matrix_, PETSC_VIEWER_STDOUT_SELF);
    //VecView(rhs_, PETSC_VIEWER_STDOUT_SELF);
    //this->view();

    matrix_changed_ = true;
    rhs_changed_ = true;
}

void LinSys_PETSC::apply_constrains( double scalar )
{
    PetscErrorCode ierr;

    // check that system matrix is assembled
    ASSERT_EQ(status_, DONE).error("System matrix and right-hand side are not assembled when applying constraints." );

    // number of constraints
    PetscInt numConstraints = static_cast<PetscInt>( constraints_.size() );

    // Additional multiplier for numerical reasons (criterion to be established)
    const PetscScalar diagScalar = static_cast<PetscScalar>( scalar );

    std::vector<PetscInt> globalDofs;
    std::vector<PetscScalar>  values;

    // Constraint iterators
    ConstraintVec_::const_iterator cIter = constraints_.begin( );
    ConstraintVec_::const_iterator cEnd  = constraints_.end( );
    // collect global dof indices and the correpsonding values
    for ( ; cIter != cEnd; ++cIter ) {
        globalDofs.push_back( static_cast<PetscInt>( cIter -> first ) );
        values.push_back( static_cast<PetscScalar>( cIter -> second ) * diagScalar );
    }

    // prepare pointers to be passed to PETSc
    PetscInt * globalDofPtr = this->makePetscPointer_( globalDofs );
    PetscScalar * valuePtr  = this->makePetscPointer_( values );

    // set matrix rows to zero 
    ierr = MatZeroRows( matrix_, numConstraints, globalDofPtr, diagScalar, PETSC_NULL, PETSC_NULL ); CHKERRV( ierr ); 
    matrix_changed_ = true;

    // set RHS entries to values (crashes if called with NULL pointers)
    if ( numConstraints ) {
        ierr = VecSetValues( rhs_, numConstraints, globalDofPtr, valuePtr, INSERT_VALUES ); CHKERRV( ierr ); 
        rhs_changed_ = true;
    }

    // perform communication in the rhs vector
    ierr = VecAssemblyBegin( rhs_ ); CHKERRV( ierr ); 
    ierr = VecAssemblyEnd(   rhs_ ); CHKERRV( ierr ); 
}


void LinSys_PETSC::set_initial_guess_nonzero(bool set_nonzero)
{
	init_guess_nonzero = set_nonzero;
}


LinSys::SolveInfo LinSys_PETSC::solve()
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
    
    chkerr(KSPCreate( comm_, &system ));
    chkerr(KSPSetOperators(system, matrix_, matrix_));


    // TODO take care of tolerances - shall we support both input file and command line petsc setting
    chkerr(KSPSetTolerances(system, r_tol_, a_tol_, PETSC_DEFAULT,PETSC_DEFAULT));
    chkerr(KSPSetTolerances(system, r_tol_, a_tol_, PETSC_DEFAULT,  max_it_));
    KSPSetFromOptions(system);
    // We set the KSP flag set_initial_guess_nonzero
    // unless KSP type is preonly.
    // In such case PETSc fails (version 3.4.1)
    if (init_guess_nonzero)
    {
    	KSPType type;
    	KSPGetType(system, &type);
    	if (strcmp(type, KSPPREONLY) != 0)
    		KSPSetInitialGuessNonzero(system, PETSC_TRUE);
    }

    {
		START_TIMER("PETSC linear solver");
		START_TIMER("PETSC linear iteration");
		chkerr(KSPSolve(system, rhs_, solution_ ));
		KSPGetConvergedReason(system,&reason);
		KSPGetIterationNumber(system,&nits);
		ADD_CALLS(nits);
    }
    // substitute by PETSc call for residual
    VecNorm(rhs_, NORM_2, &residual_norm_);
    
    LogOut().fmt("convergence reason {}, number of iterations is {}\n", reason, nits);

    // get residual norm
    KSPGetResidualNorm(system, &solution_precision_);

    // TODO: I do not understand this 
    //Profiler::instance()->set_timer_subframes("SOLVING MH SYSTEM", nits);

    chkerr(KSPDestroy(&system));

    return LinSys::SolveInfo(static_cast<int>(reason), static_cast<int>(nits));

}

void LinSys_PETSC::view(string text )
{
    FilePath matFileName(text + "_flow123d_matrix.m",FilePath::FileType::output_file);
    FilePath rhsFileName(text + "_flow123d_rhs.m",FilePath::FileType::output_file);
    FilePath solFileName(text + "_flow123d_sol.m",FilePath::FileType::output_file);

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
}

LinSys_PETSC::~LinSys_PETSC( )
{
    if (matrix_ != NULL) { chkerr(MatDestroy(&matrix_)); }
    chkerr(VecDestroy(&rhs_));

    if (residual_ != NULL) chkerr(VecDestroy(&residual_));
    if (v_rhs_ != NULL) delete[] v_rhs_;
}



void LinSys_PETSC::set_from_input(const Input::Record in_rec)
{
	LinSys::set_from_input( in_rec );

	// PETSC specific parameters
    // If parameters are specified in input file, they are used,
    // otherwise keep settings provided in constructor of LinSys_PETSC.
    std::string user_params = in_rec.val<string>("options");
	if (user_params != "") params_ = user_params;
}


double LinSys_PETSC::get_solution_precision()
{
	return solution_precision_;
}


double LinSys_PETSC::compute_residual()
{
    MatMult(matrix_, solution_, residual_);
    VecAXPY(residual_,-1.0, rhs_);
    double residual_norm;
    VecNorm(residual_, NORM_2, &residual_norm);
    return residual_norm;
}
