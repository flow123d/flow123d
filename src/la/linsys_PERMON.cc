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


//#include <boost/bind.hpp>

namespace it = Input::Type;

const it::Record & LinSys_PERMON::get_input_type() {
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


const int LinSys_PERMON::registrar = LinSys_PERMON::get_input_type().size();


LinSys_PERMON::LinSys_PERMON( const Distribution * rows_ds, const std::string &params)
        : LinSys_PETSC( rows_ds ),
          params_(params)
{
    matrix_ineq_ = NULL;
    ineq_ = NULL;
}

LinSys_PERMON::LinSys_PERMON( LinSys_PERMON &other )
	: LinSys_PETSC(other), params_(other.params_)
{
	MatCopy(other.matrix_ineq_, matrix_ineq_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.ineq_, ineq_);
}

void LinSys_PERMON::set_inequality(Mat matrix_ineq, Vec ineq)
{
  // TODO ref count?
  matrix_ineq_ = matrix_ineq;
  ineq_ = ineq;
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
    
    chkerr(QPCreate(comm_, &system));
    chkerr(QPSetOperator(system, matrix_));
    chkerr(QPSetRhs(system, rhs_));
    chkerr(QPSetInitialVector(system, solution_));
    if (ineq_) {
      //chkerr(MatScale(matrix_ineq_,-1.));
      //chkerr(VecScale(ineq_,-1.));
      // Bx<=c
      chkerr(QPSetIneq(system, matrix_ineq_, ineq_));
      chkerr(QPSetIneq(system, matrix_ineq_, ineq_));
      chkerr(QPTDualize(system, MAT_INV_MONOLITHIC, MAT_REG_NONE));
    }
    // Set runtime options, e.g -qp_chain_view_kkt
    chkerr(QPSetFromOptions(system));
      
    
    chkerr(QPSCreate(comm_, &solver));
    chkerr(QPSSetQP(solver, system));

    // TODO take care of tolerances - shall we support both input file and command line petsc setting
    chkerr(QPSSetTolerances(solver, r_tol_, a_tol_, PETSC_DEFAULT,PETSC_DEFAULT));
    chkerr(QPSSetTolerances(solver, r_tol_, a_tol_, PETSC_DEFAULT,  max_it_));
    chkerr(QPSSetFromOptions(solver));

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
    //KSPGetResidualNorm(system, &solution_precision_);

    // TODO: I do not understand this 
    //Profiler::instance()->set_timer_subframes("SOLVING MH SYSTEM", nits);

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
    // TODO cleanup ineq objects
}



void LinSys_PERMON::set_from_input(const Input::Record in_rec)
{
	LinSys::set_from_input( in_rec );

	// PETSC specific parameters
    // If parameters are specified in input file, they are used,
    // otherwise keep settings provided in constructor of LinSys_PERMON.
    std::string user_params = in_rec.val<string>("options");
	if (user_params != "") params_ = user_params;
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
    double residual_norm;
    VecNorm(residual_, NORM_2, &residual_norm);
    return residual_norm;
}
