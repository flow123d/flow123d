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
#include "permonmat.h"
#include <cmath>


//#include <boost/bind.hpp>

namespace it = Input::Type;

namespace {

/// Solve a configured QP using the common PERMON tolerances and iteration accounting.
LinSys::SolveInfo solve_permon_qp(QP qp, QPS *qps, MPI_Comm comm, double r_tol, double a_tol,
                                  double d_tol, int max_it, const char *prefix,
                                  const char *type)
{
    chkerr(QPSCreate(comm, qps));
    if (prefix) chkerr(QPSSetOptionsPrefix(*qps, prefix));
    chkerr(QPSSetQP(*qps, qp));
    if (type) chkerr(QPSSetType(*qps, type));
    chkerr(QPSSetTolerances(*qps, r_tol, a_tol, d_tol, max_it));
    chkerr(QPSSetFromOptions(*qps));
    KSPConvergedReason reason;
    PetscInt nits;
    {
        START_TIMER("PERMON linear solver");
        START_TIMER("PERMON linear iteration");
        chkerr(QPSSolve(*qps));
        chkerr(QPSGetConvergedReason(*qps, &reason));
        chkerr(QPSGetIterationNumber(*qps, &nits));
        ADD_CALLS(nits);
    }
    return LinSys::SolveInfo(static_cast<int>(reason), static_cast<int>(nits));
}

// Own temporary PETSc/PERMON objects through every return path. Input objects
// and the LinSys_PERMON multiplier fields remain outside this ownership scope.
struct DualWork {
    Mat nest = nullptr, constraints = nullptr, constraints_transpose = nullptr;
    Mat inverse = nullptr, hessian = nullptr;
    Vec multipliers = nullptr, offsets = nullptr, rhs = nullptr;
    Vec lower = nullptr, upper = nullptr;
    Vec primal_work = nullptr, shifted_rhs = nullptr;
    Mat active_box = nullptr;
    IS active_rows = nullptr, columns = nullptr;
    Vec active_lower = nullptr, active_upper = nullptr, active_offset = nullptr;
    QP qp = nullptr;
    QPS qps = nullptr;

    ~DualWork() {
        QPSDestroy(&qps); QPDestroy(&qp);
        VecDestroy(&multipliers); VecDestroy(&offsets); VecDestroy(&rhs);
        VecDestroy(&lower); VecDestroy(&upper); VecDestroy(&primal_work); VecDestroy(&shifted_rhs);
        MatDestroy(&hessian); MatDestroy(&constraints_transpose); MatDestroy(&inverse);
        MatDestroy(&constraints); MatDestroy(&nest);
        MatDestroy(&active_box); ISDestroy(&active_rows); ISDestroy(&columns);
        VecDestroy(&active_lower); VecDestroy(&active_upper); VecDestroy(&active_offset);
    }
};

// Return an owned copy; VecGetSubVector may otherwise expose a borrowed view.
void copy_selected_rows(Vec source, IS rows, Vec *selected)
{
    Vec view;
    chkerr(VecGetSubVector(source, rows, &view));
    chkerr(VecDuplicate(view, selected));
    chkerr(VecCopy(view, *selected));
    chkerr(VecRestoreSubVector(source, rows, &view));
}

} // namespace

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
        : LinSys_PETSC( rows_ds, params )
{
    matrix_ineq_ = NULL;
    ineq_ = NULL;
}

void LinSys_PERMON::set_inequality(Mat matrix_ineq, Vec ineq)
{
  // TODO ref count?
  matrix_ineq_ = matrix_ineq;
  ineq_ = ineq;
  chkerr(VecDestroy(&ineq_multipliers_));
}


void LinSys_PERMON::set_dual_box(Mat B, Vec offset, Vec lower, Vec upper)
{
    PetscInt rows, columns, lower_size, upper_size, offset_size;
    chkerr(MatGetLocalSize(B, &rows, &columns));
    chkerr(VecGetLocalSize(lower, &lower_size));
    chkerr(VecGetLocalSize(upper, &upper_size));
    chkerr(VecGetLocalSize(offset, &offset_size));
    ASSERT_PERMANENT_EQ(rows, lower_size).error("Dual lower bounds must use the B row layout.");
    ASSERT_PERMANENT_EQ(rows, upper_size).error("Dual upper bounds must use the B row layout.");
    ASSERT_PERMANENT_EQ(rows, offset_size).error("Dual offset must use the B row layout.");
    ASSERT_PERMANENT_EQ(columns, static_cast<PetscInt>(rows_ds_->lsize()))
        .error("B must use the original primal column layout.");
    chkerr(PetscObjectReference((PetscObject)B));
    chkerr(PetscObjectReference((PetscObject)offset));
    chkerr(PetscObjectReference((PetscObject)lower));
    chkerr(PetscObjectReference((PetscObject)upper));
    chkerr(MatDestroy(&matrix_dual_box_));
    chkerr(VecDestroy(&dual_box_lower_));
    chkerr(VecDestroy(&dual_box_upper_));
    chkerr(VecDestroy(&dual_box_offset_));
    matrix_dual_box_ = B;
    dual_box_lower_ = lower;
    dual_box_upper_ = upper;
    dual_box_offset_ = offset;
    chkerr(VecDestroy(&ineq_multipliers_));
    chkerr(VecDestroy(&dual_box_multipliers_));
}


LinSys::SolveInfo LinSys_PERMON::solve_dual_qp()
{
    // Allocate temporary storage for the dual of the original primal problem:
    // min 0.5*x^T*K*x - f^T*x +
    // max_{lower <= q <= upper} q^T*(B*x-offset),
    // subject to A*x <= b. K must be positive definite. The multiplier order
    // below is [lambda for A; active q for B], with lambda >= 0 and box-bounded q.
    // All B columns and recovered x use the original global primal numbering.
    DualWork dual;

    // Read the locally owned box bounds and offsets for validation and reduction.
    PetscInt local_rows;
    const PetscScalar *lower_bounds, *upper_bounds, *offset;
    chkerr(VecGetLocalSize(dual_box_upper_, &local_rows));
    chkerr(VecGetArrayRead(dual_box_lower_, &lower_bounds));
    chkerr(VecGetArrayRead(dual_box_upper_, &upper_bounds));
    chkerr(VecGetArrayRead(dual_box_offset_, &offset));

    // Reset the stored multipliers and allocate the full B-row result vector.
    chkerr(VecDestroy(&ineq_multipliers_));
    chkerr(VecDestroy(&dual_box_multipliers_));
    chkerr(VecDuplicate(dual_box_upper_, &dual_box_multipliers_));
    PetscScalar *fixed_multipliers;
    chkerr(VecGetArray(dual_box_multipliers_, &fixed_multipliers));

    // Validate each local row, retain fixed multipliers, and collect global B row
    // indices with lower < upper for the reduced dual problem.
    int local_error = 0, global_error;
    PetscInt first_box_row;
    chkerr(MatGetOwnershipRange(matrix_dual_box_, &first_box_row, nullptr));
    std::vector<PetscInt> active_rows;
    for (PetscInt i = 0; i < local_rows; ++i) {
        const PetscScalar lower = lower_bounds[i], upper = upper_bounds[i];
        if (std::isnan(lower) || std::isnan(upper) || lower > upper
            || lower == PETSC_INFINITY || upper == PETSC_NINFINITY
            || (!std::isfinite(lower) && lower > 0.0)
            || (!std::isfinite(upper) && upper < 0.0)) local_error |= 1;
        if (!std::isfinite(offset[i])) local_error |= 2;
        fixed_multipliers[i] = lower == upper ? lower : 0.0;
        if (lower < upper) active_rows.push_back(first_box_row + i);
    }

    // Restore the borrowed arrays before checking validation errors collectively.
    chkerr(VecRestoreArray(dual_box_multipliers_, &fixed_multipliers));
    chkerr(VecRestoreArrayRead(dual_box_lower_, &lower_bounds));
    chkerr(VecRestoreArrayRead(dual_box_upper_, &upper_bounds));
    chkerr(VecRestoreArrayRead(dual_box_offset_, &offset));

    // Reject invalid bounds or nonfinite offsets on every MPI rank.
    MPI_Allreduce(&local_error, &global_error, 1, MPI_INT, MPI_BOR, comm_);
    ASSERT_PERMANENT((global_error & 1) == 0)
        .error("Dual box bounds must satisfy lower <= upper and admit finite multipliers.");
    ASSERT_PERMANENT((global_error & 2) == 0).error("Dual box offset must be finite.");

    // Eliminate fixed multipliers by shifting f to f - B_fixed^T*q_fixed;
    // avoid passing coincident bounds to MPGP, which handles them unreliably.
    chkerr(VecDuplicate(rhs_, &dual.shifted_rhs));
    chkerr(MatMultTranspose(matrix_dual_box_, dual_box_multipliers_, dual.shifted_rhs));
    chkerr(VecAYPX(dual.shifted_rhs, -1.0, rhs_));

    // Build the distributed selection of nonfixed B rows and determine which
    // multiplier blocks remain in the reduced dual problem.
    chkerr(ISCreateGeneral(comm_, active_rows.size(), active_rows.data(), PETSC_COPY_VALUES, &dual.active_rows));
    PetscInt n_active, first_column, last_column, n_ineq = 0;
    chkerr(ISGetSize(dual.active_rows, &n_active));
    if (matrix_ineq_) chkerr(MatGetSize(matrix_ineq_, &n_ineq, nullptr));
    const bool has_box = n_active > 0, has_ineq = n_ineq > 0;

    // Extract the nonfixed box rows while retaining all original primal columns.
    if (has_box) {
        // Select original global primal columns explicitly for MatCreateSubMatrix.
        chkerr(MatGetOwnershipRangeColumn(matrix_dual_box_, &first_column, &last_column));
        chkerr(ISCreateStride(comm_, last_column-first_column, first_column, 1, &dual.columns));
        chkerr(MatCreateSubMatrix(matrix_dual_box_, dual.active_rows, dual.columns,
                                 MAT_INITIAL_MATRIX, &dual.active_box));

        // Copy the corresponding bounds and offsets into owned reduced vectors.
        copy_selected_rows(dual_box_lower_, dual.active_rows, &dual.active_lower);
        copy_selected_rows(dual_box_upper_, dual.active_rows, &dual.active_upper);
        copy_selected_rows(dual_box_offset_, dual.active_rows, &dual.active_offset);
    }

    // Configure the monolithic K^-1 operator without regularization, using the
    // positive-definite original primal Hessian and dual_inverse_ options.
    chkerr(MatCreateInv(matrix_, MAT_INV_MONOLITHIC, &dual.inverse));
    chkerr(MatSetOptionsPrefix(dual.inverse, "dual_inverse_"));
    chkerr(MatInvSetRegularizationType(dual.inverse, MAT_REG_NONE));
    chkerr(MatSetFromOptions(dual.inverse));
    chkerr(MatInvSetUp(dual.inverse));

    // Solve the shifted primal directly when no nonfixed multipliers remain.
    if (!has_ineq && !has_box) {
        chkerr(MatMult(dual.inverse, dual.shifted_rhs, solution_));
        solution_precision_ = compute_residual();
        return LinSys::SolveInfo(KSP_CONVERGED_ATOL, 0);
    }

    // Assemble D = [A; B_active], omitting absent blocks, and borrow the MatNest
    // row index sets for packing and recovery. Keep the nest alive: converted
    // AIJ rows need not have a simple global block offset.
    const PetscInt n_blocks = static_cast<PetscInt>(has_ineq) + static_cast<PetscInt>(has_box);
    Mat blocks[2] = {has_ineq ? matrix_ineq_ : dual.active_box, dual.active_box};
    IS block_rows[2] = {nullptr, nullptr};
    chkerr(MatCreateNest(comm_, n_blocks, nullptr, 1, nullptr, blocks, &dual.nest));
    chkerr(MatNestGetISs(dual.nest, block_rows, nullptr));
    chkerr(MatConvert(dual.nest, MATAIJ, MAT_INITIAL_MATRIX, &dual.constraints));

    // Allocate primal work and dual vectors in the column and row layouts of D.
    chkerr(MatCreateVecs(dual.constraints, &dual.primal_work, &dual.multipliers));
    chkerr(VecDuplicate(dual.multipliers, &dual.offsets));
    chkerr(VecDuplicate(dual.multipliers, &dual.rhs));
    chkerr(VecDuplicate(dual.multipliers, &dual.lower));
    chkerr(VecDuplicate(dual.multipliers, &dual.upper));

    // Start from zero multipliers and unbounded intervals before filling blocks.
    chkerr(VecSet(dual.multipliers, 0.0));
    chkerr(VecSet(dual.lower, PETSC_NINFINITY));
    chkerr(VecSet(dual.upper, PETSC_INFINITY));

    // Pack [b; offset_active] and assign lambda >= 0 or the supplied box bounds.
    for (PetscInt block = 0; block < n_blocks; ++block) {
        const bool is_box_block = has_box && block == n_blocks - 1;
        Vec block_offset, lower, upper;
        chkerr(VecGetSubVector(dual.offsets, block_rows[block], &block_offset));
        chkerr(VecGetSubVector(dual.lower, block_rows[block], &lower));
        chkerr(VecGetSubVector(dual.upper, block_rows[block], &upper));
        chkerr(VecCopy(is_box_block ? dual.active_offset : ineq_, block_offset));
        if (is_box_block) {
            chkerr(VecCopy(dual.active_lower, lower));
            chkerr(VecCopy(dual.active_upper, upper));
        } else {
            chkerr(VecSet(lower, 0.0));
        }
        chkerr(VecRestoreSubVector(dual.offsets, block_rows[block], &block_offset));
        chkerr(VecRestoreSubVector(dual.lower, block_rows[block], &lower));
        chkerr(VecRestoreSubVector(dual.upper, block_rows[block], &upper));
    }

    // Construct the dual Hessian D*K^-1*D^T as an implicit matrix product.
    chkerr(MatTranspose(dual.constraints, MAT_INITIAL_MATRIX, &dual.constraints_transpose));
    // List factors in reverse multiplication order for PERMON MatCreateProd.
    Mat factors[3] = {dual.constraints_transpose, dual.inverse, dual.constraints};
    chkerr(MatCreateProd(comm_, 3, factors, &dual.hessian));

    // Form the dual RHS D*K^-1*f_shifted - offsets.
    chkerr(MatMult(dual.inverse, dual.shifted_rhs, dual.primal_work));
    chkerr(MatMult(dual.constraints, dual.primal_work, dual.rhs));
    chkerr(VecAXPY(dual.rhs, -1.0, dual.offsets));

    // Define min_z 0.5*z^T*hessian*z - rhs^T*z with the assembled multiplier bounds.
    chkerr(QPCreate(comm_, &dual.qp));
    chkerr(QPSetOptionsPrefix(dual.qp, "dual_"));
    chkerr(QPSetOperator(dual.qp, dual.hessian));
    chkerr(QPSetRhs(dual.qp, dual.rhs));
    chkerr(QPSetInitialVector(dual.qp, dual.multipliers));
    chkerr(QPSetBox(dual.qp, nullptr, dual.lower, dual.upper));

    // Solve with MPGP defaults and allow dual_ runtime options to override them.
    LinSys::SolveInfo info = solve_permon_qp(dual.qp, &dual.qps, comm_, r_tol_, a_tol_, d_tol_,
                                              max_it_, "dual_", QPSMPGP);

    // Borrow z from the unscaled dual after PERMON post-solve has mapped any
    // transformed child solution back to this QP's distributed row numbering.
    Vec multipliers;
    chkerr(QPGetSolutionVector(dual.qp, &multipliers));

    // Recover x = K^-1*(f_shifted - D^T*z) in original global primal numbering.
    chkerr(MatMultTranspose(dual.constraints, multipliers, dual.primal_work));
    chkerr(VecAYPX(dual.primal_work, -1.0, dual.shifted_rhs));
    chkerr(MatMult(dual.inverse, dual.primal_work, solution_));

    // Scatter the solved blocks into the full A-row and B-row multiplier vectors;
    // preserve the fixed B-row values stored before the reduction.
    for (PetscInt block = 0; block < n_blocks; ++block) {
        const bool is_box_block = has_box && block == n_blocks - 1;
        Vec sub;
        chkerr(VecGetSubVector(multipliers, block_rows[block], &sub));
        if (is_box_block) {
            Vec active_multipliers;
            chkerr(VecGetSubVector(dual_box_multipliers_, dual.active_rows, &active_multipliers));
            chkerr(VecCopy(sub, active_multipliers));
            chkerr(VecRestoreSubVector(dual_box_multipliers_, dual.active_rows, &active_multipliers));
        } else {
            chkerr(VecDuplicate(ineq_, &ineq_multipliers_));
            chkerr(VecCopy(sub, ineq_multipliers_));
        }
        chkerr(VecRestoreSubVector(multipliers, block_rows[block], &sub));
    }

    // Report the original primal's algebraic KKT norm, combining equilibrium,
    // inequality complementarity, and box projection residuals.
    solution_precision_ = compute_residual();
    LogOut().fmt("Dual QP: reason {}, iterations {}, algebraic KKT residual {}\n",
                 info.converged_reason, info.n_iterations, solution_precision_);
    return info;
}

LinSys::SolveInfo LinSys_PERMON::solve()
{
    // Discard multipliers from the previous solve, including now-empty constraints.
    chkerr(VecDestroy(&ineq_multipliers_));
    chkerr(VecDestroy(&dual_box_multipliers_));
    const char *petsc_dflt_opt;
    
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
    chkerr(PetscOptionsInsertString(NULL, params_.c_str())); // overwrites previous options values
    
    chkerr(MatSetOption( matrix_, MAT_USE_INODES, PETSC_FALSE ));

    if (matrix_dual_box_) {
        PetscInt n_box;
        chkerr(MatGetSize(matrix_dual_box_, &n_box, nullptr));
        if (n_box) return solve_dual_qp();
    }
    
    chkerr(QPCreate(comm_, &system));
    chkerr(QPSetOperator(system, matrix_));
    chkerr(QPSetRhs(system, rhs_));
    chkerr(QPSetInitialVector(system, solution_));
    QP unscaled_dual = nullptr;
    PetscInt n_ineq = 0;
    if (matrix_ineq_ && ineq_) chkerr(MatGetSize(matrix_ineq_, &n_ineq, nullptr));
    if (n_ineq > 0) {
        // Keep the direct dual of the original primal A*x <= b. Its solution
        // uses original A rows after post-solve, even with transformed children.
        chkerr(QPSetIneq(system, matrix_ineq_, ineq_));
        chkerr(QPTDualize(system, MAT_INV_MONOLITHIC, MAT_REG_NONE));
        chkerr(QPGetChild(system, &unscaled_dual));
    }
    // Set runtime options, e.g -qp_chain_view_kkt
    chkerr(QPSetFromOptions(system));
      
    
    LinSys::SolveInfo info = solve_permon_qp(system, &solver, comm_, r_tol_, a_tol_, d_tol_,
                                              max_it_, nullptr, nullptr);
    if (unscaled_dual) {
        // Copy the borrowed physical multipliers before destroying the QP chain.
        Vec multipliers;
        chkerr(QPGetSolutionVector(unscaled_dual, &multipliers));
        PetscInt global_size, local_size, expected_local_size;
        chkerr(VecGetSize(multipliers, &global_size));
        chkerr(VecGetLocalSize(multipliers, &local_size));
        chkerr(VecGetLocalSize(ineq_, &expected_local_size));
        ASSERT_PERMANENT_EQ(global_size, n_ineq).error("Contact multipliers must use original inequality rows.");
        ASSERT_PERMANENT_EQ(local_size, expected_local_size).error("Contact multiplier ownership must match A rows.");
        chkerr(VecDuplicate(ineq_, &ineq_multipliers_));
        chkerr(VecCopy(multipliers, ineq_multipliers_));
    }
    // substitute by PETSc call for residual
    //VecNorm(rhs_, NORM_2, &residual_norm_);
    
    LogOut().fmt("convergence reason {}, number of iterations is {}\n",
                 info.converged_reason, info.n_iterations);

    // get residual norm
    //KSPGetResidualNorm(system, &solution_precision_);

    // TODO: I do not understand this 
    //Profiler::instance()->set_timer_subframes("SOLVING MH SYSTEM", nits);

    chkerr(QPSDestroy(&solver));
    chkerr(QPDestroy(&system));

    return info;

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
    MatDestroy(&matrix_dual_box_);
    VecDestroy(&dual_box_lower_);
    VecDestroy(&dual_box_upper_);
    VecDestroy(&dual_box_offset_);
    VecDestroy(&ineq_multipliers_);
    VecDestroy(&dual_box_multipliers_);
}




double LinSys_PERMON::get_solution_precision()
{
	return solution_precision_;
}


double LinSys_PERMON::compute_residual()
{
    if (dual_box_multipliers_) {
        struct Work {
            Vec column = nullptr, inequality = nullptr, zero = nullptr, box = nullptr;
            ~Work() {
                VecDestroy(&column); VecDestroy(&inequality); VecDestroy(&zero); VecDestroy(&box);
            }
        } w;
        chkerr(VecDuplicate(solution_, &w.column));
        chkerr(MatMult(matrix_, solution_, residual_));
        chkerr(VecAXPY(residual_, -1.0, rhs_));
        if (ineq_multipliers_) {
            chkerr(MatMultTranspose(matrix_ineq_, ineq_multipliers_, w.column));
            chkerr(VecAXPY(residual_, 1.0, w.column));
        }
        chkerr(MatMultTranspose(matrix_dual_box_, dual_box_multipliers_, w.column));
        chkerr(VecAXPY(residual_, 1.0, w.column));
        PetscReal equilibrium, inequality = 0.0, box;
        chkerr(VecNorm(residual_, NORM_2, &equilibrium));
        if (ineq_multipliers_) {
            chkerr(VecDuplicate(ineq_, &w.inequality));
            chkerr(VecDuplicate(ineq_, &w.zero));
            chkerr(VecSet(w.zero, 0.0));
            chkerr(MatMult(matrix_ineq_, solution_, w.inequality));
            chkerr(VecAXPY(w.inequality, -1.0, ineq_));
            chkerr(VecAXPY(w.inequality, 1.0, ineq_multipliers_));
            chkerr(VecPointwiseMax(w.inequality, w.inequality, w.zero));
            chkerr(VecAYPX(w.inequality, -1.0, ineq_multipliers_));
            chkerr(VecNorm(w.inequality, NORM_2, &inequality));
        }
        chkerr(VecDuplicate(dual_box_upper_, &w.box));
        chkerr(MatMult(matrix_dual_box_, solution_, w.box));
        chkerr(VecAXPY(w.box, -1.0, dual_box_offset_));
        chkerr(VecAXPY(w.box, 1.0, dual_box_multipliers_));
        chkerr(VecPointwiseMin(w.box, w.box, dual_box_upper_));
        chkerr(VecPointwiseMax(w.box, w.box, dual_box_lower_));
        chkerr(VecAYPX(w.box, -1.0, dual_box_multipliers_));
        chkerr(VecNorm(w.box, NORM_2, &box));
        return std::hypot(equilibrium, std::hypot(inequality, box));
    }
    // TODO return ||A*x - b + B'*lambda||?
    MatMult(matrix_, solution_, residual_);
    VecAXPY(residual_,-1.0, rhs_);
    double residual_norm;
    VecNorm(residual_, NORM_2, &residual_norm);
    return residual_norm;
}
