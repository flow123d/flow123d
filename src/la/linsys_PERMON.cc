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
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>
#include <vector>
#include "petscvec.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscblaslapack.h"
#include "permonvec.h"
#include "permonmat.h"
#include "system/sys_profiler.hh"
#include "system/system.hh"
#include "fem/dofhandler.hh"







namespace it = Input::Type;

namespace {

// Auxiliary function for converting matrix from undecomposed to AIJ indexing.
PetscErrorCode convert_mat_is_to_aij(Mat matrix_is, Mat *matrix_aij)
{
    PetscFunctionBegin;
    PetscCheck(matrix_is, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null inequality matrix.");
    PetscCheck(matrix_aij, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null output inequality matrix pointer.");
    PetscCall(MatConvert(matrix_is, MATAIJ, MAT_INITIAL_MATRIX, matrix_aij));
    PetscFunctionReturn(PETSC_SUCCESS);
}

// Auxiliary function for converting IS from undecomposed to decomposed indexing.
// input:
//   matrix_is: Matrix to be converted.
//   isnz: Index set of nonzero rows/columns in Hessian.
//   A_decomposed: Hessian from which the column distribution is taken.
// output:
//   matrix_decomp: Resulting matrix in decomposed indexing.
PetscErrorCode convert_mat_is_to_feti_decomposed(Mat matrix_is,
                                                  IS isnz,
                                                  Mat A_decomposed,
                                                  Mat *matrix_decomp)
{
    Mat Bloc = NULL, Bloc_reduced = NULL;
    IS ris_all = NULL;
    ISLocalToGlobalMapping row_l2g = NULL, col_l2g_dec = NULL;
    const PetscInt *row_gidx = NULL, *col_gidx_dec = NULL;
    PetscInt nrows_local, ncols_local_dec;
    PetscInt M_rows, N_cols_dec;

    PetscFunctionBegin;
    PetscCheck(matrix_is, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null inequality matrix.");
    PetscCheck(isnz, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null reduced-column index set.");
    PetscCheck(A_decomposed, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null decomposed operator.");
    PetscCheck(matrix_decomp, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null output inequality matrix pointer.");

    PetscCall(MatISGetLocalMat(matrix_is, &Bloc));
    PetscCall(MatGetLocalSize(Bloc, &nrows_local, NULL));
    PetscCall(ISCreateStride(PETSC_COMM_SELF, nrows_local, 0, 1, &ris_all));
    PetscCall(MatCreateSubMatrix(Bloc, ris_all, isnz, MAT_INITIAL_MATRIX, &Bloc_reduced));

    PetscCall(MatGetLocalSize(Bloc_reduced, &nrows_local, &ncols_local_dec));
    PetscCall(MatGetSize(matrix_is, &M_rows, NULL));
    PetscCall(MatGetSize(A_decomposed, NULL, &N_cols_dec));

    PetscCall(MatGetLocalToGlobalMapping(matrix_is, &row_l2g, NULL));
    PetscCall(ISLocalToGlobalMappingGetIndices(row_l2g, &row_gidx));
    PetscCall(MatGetLocalToGlobalMapping(A_decomposed, &col_l2g_dec, NULL));
    PetscCall(ISLocalToGlobalMappingGetIndices(col_l2g_dec, &col_gidx_dec));

    PetscCall(MatCreateAIJ(PetscObjectComm((PetscObject)matrix_is),
                           nrows_local, ncols_local_dec,
                           M_rows,     N_cols_dec,
                           32, NULL, 32, NULL,
                           matrix_decomp));

    for (PetscInt i = 0; i < nrows_local; ++i) {
        PetscInt ncols;
        const PetscInt *cols_loc;
        const PetscScalar *vals;
        PetscInt row_global = row_gidx[i];

        PetscCall(MatGetRow(Bloc_reduced, i, &ncols, &cols_loc, &vals));
        if (ncols > 0) {
            std::vector<PetscInt> cols_global(ncols);
            for (PetscInt j = 0; j < ncols; ++j) {
                PetscInt loc_col = cols_loc[j];
                PetscCheck(loc_col >= 0 && loc_col < ncols_local_dec,
                           PetscObjectComm((PetscObject)matrix_is), PETSC_ERR_ARG_OUTOFRANGE,
                           "Reduced local inequality column index out of range.");
                cols_global[j] = col_gidx_dec[loc_col];
            }
            PetscCall(MatSetValues(*matrix_decomp, 1, &row_global, ncols, cols_global.data(), vals, ADD_VALUES));
        }
        PetscCall(MatRestoreRow(Bloc_reduced, i, &ncols, &cols_loc, &vals));
    }

    PetscCall(MatAssemblyBegin(*matrix_decomp, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*matrix_decomp, MAT_FINAL_ASSEMBLY));
    PetscCall(ISLocalToGlobalMappingRestoreIndices(col_l2g_dec, &col_gidx_dec));
    PetscCall(ISLocalToGlobalMappingRestoreIndices(row_l2g, &row_gidx));
    PetscCall(ISDestroy(&ris_all));
    PetscCall(MatDestroy(&Bloc_reduced));
    PetscCall(MatISRestoreLocalMat(matrix_is, &Bloc));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode create_reduced_nullspace_local(Mat nullspace_is,
                                               IS isnz,
                                               Mat *Rloc_reduced)
{
    Mat Rloc = NULL;
    IS cis_all = NULL;
    PetscInt nrows_local, ncols_local;

    PetscFunctionBegin;
    PetscCheck(nullspace_is, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Null operator nullspace matrix.");
    PetscCheck(isnz, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Null reduced-row index set.");
    PetscCheck(Rloc_reduced, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Null output reduced local nullspace matrix pointer.");

    PetscCall(MatISGetLocalMat(nullspace_is, &Rloc));
    PetscCall(MatGetLocalSize(Rloc, &nrows_local, &ncols_local));
    PetscCall(ISCreateStride(PETSC_COMM_SELF, ncols_local, 0, 1, &cis_all));
    PetscCall(MatCreateSubMatrix(Rloc, isnz, cis_all,
                                 MAT_INITIAL_MATRIX, Rloc_reduced));
    PetscCall(MatISRestoreLocalMat(nullspace_is, &Rloc));
    PetscCall(ISDestroy(&cis_all));

    PetscCall(MatAssemblyBegin(*Rloc_reduced, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*Rloc_reduced, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode convert_nullspace_to_feti_blockdiag(Mat Rloc_reduced,
                                                   Mat *nullspace_blockdiag)
{
    PetscFunctionBegin;
    PetscCheck(Rloc_reduced, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Null reduced local nullspace matrix.");
    PetscCheck(nullspace_blockdiag, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Null output nullspace matrix pointer.");
    PetscCall(MatCreateBlockDiag(PETSC_COMM_WORLD,
                                 Rloc_reduced,
                                 nullspace_blockdiag));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode print_feti_user_nullspace_dimension(MPI_Comm comm, Mat R)
{
    PetscBool enabled = PETSC_FALSE;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_user_nullspace_dimension",
                                  &enabled, NULL));
    if (!enabled) PetscFunctionReturn(PETSC_SUCCESS);

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    PetscInt m_local, n_local, m_global, n_global;
    PetscCall(MatGetLocalSize(R, &m_local, &n_local));
    PetscCall(MatGetSize(R, &m_global, &n_global));

    PetscCall(PetscSynchronizedPrintf(comm,
        "[feti user nullspace] rank: %d, R_local_size: %d x %d, "
        "R_global_size: %d x %d, local_dim: %d, global_dim: %d\n",
        rank, (int)m_local, (int)n_local, (int)m_global, (int)n_global,
        (int)n_local, (int)n_global));
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode print_feti_user_nullspace_test(MPI_Comm comm, Mat Aloc_feti, Mat Rloc_reduced)
{
    PetscBool enabled = PETSC_FALSE, dimension_enabled = PETSC_FALSE;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_user_nullspace_test",
                                  &enabled, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_user_nullspace_dimension",
                                  &dimension_enabled, NULL));
    if (!enabled && !dimension_enabled) PetscFunctionReturn(PETSC_SUCCESS);

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    PetscInt Am, An, Rm, Rn;
    PetscCall(MatGetSize(Aloc_feti, &Am, &An));
    PetscCall(MatGetSize(Rloc_reduced, &Rm, &Rn));
    PetscCheck(Am == An && An == Rm, comm, PETSC_ERR_ARG_SIZ,
               "Incompatible local FETI Hessian and nullspace sizes.");

    Mat AR = NULL;
    PetscCall(MatMatMult(Aloc_feti, Rloc_reduced, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AR));

    PetscReal norm_A = 0.0, norm_R = 0.0, norm_AR = 0.0;
    PetscCall(MatNorm(Aloc_feti, NORM_FROBENIUS, &norm_A));
    PetscCall(MatNorm(Rloc_reduced, NORM_FROBENIUS, &norm_R));
    PetscCall(MatNorm(AR, NORM_FROBENIUS, &norm_AR));
    PetscCall(MatDestroy(&AR));

    PetscReal residual_rel = norm_AR;
    if (norm_A > 0.0 && norm_R > 0.0) residual_rel /= norm_A * norm_R;

    Mat Gram = NULL, dense = NULL;
    PetscCall(MatTransposeMatMult(Rloc_reduced, Rloc_reduced, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Gram));
    PetscCall(MatConvert(Gram, MATSEQDENSE, MAT_INITIAL_MATRIX, &dense));
    PetscCall(MatDestroy(&Gram));

    PetscScalar *array = NULL;
    PetscCall(MatDenseGetArray(dense, &array));

    PetscBLASInt bn, lda, info;
    PetscCall(PetscBLASIntCast(Rn, &bn));
    lda = bn;

    PetscReal *eigenvalues = NULL;
    PetscCall(PetscMalloc1(Rn, &eigenvalues));

    PetscBLASInt lwork = -1;
    PetscScalar work_query;
    LAPACKsyev_("N", "U", &bn, array, &lda, eigenvalues,
                &work_query, &lwork, &info);
    PetscCheck(info == 0, comm, PETSC_ERR_LIB,
               "LAPACKsyev workspace query failed with info=%d", (int)info);

    lwork = static_cast<PetscBLASInt>(PetscRealPart(work_query));
    PetscScalar *work = NULL;
    PetscCall(PetscMalloc1(lwork, &work));

    LAPACKsyev_("N", "U", &bn, array, &lda, eigenvalues,
                work, &lwork, &info);

    PetscCall(MatDenseRestoreArray(dense, &array));
    PetscCall(MatDestroy(&dense));
    PetscCall(PetscFree(work));

    PetscCheck(info == 0, comm, PETSC_ERR_LIB,
               "LAPACKsyev failed with info=%d", (int)info);

    PetscReal lambda_max = 0.0;
    for (PetscInt i = 0; i < Rn; ++i)
        lambda_max = std::max(lambda_max, PetscAbsReal(eigenvalues[i]));

    PetscReal rank_rel_tol = 1e-10;
    PetscCall(PetscOptionsGetReal(NULL, NULL,
                                  "-flow123d_feti_user_nullspace_rank_rel_tol",
                                  &rank_rel_tol, NULL));
    PetscReal rank_tol = rank_rel_tol * (lambda_max > 0.0 ? lambda_max : 1.0);
    PetscInt rank_estimate = 0;
    for (PetscInt i = 0; i < Rn; ++i)
        if (eigenvalues[i] > rank_tol) ++rank_estimate;

    PetscCall(PetscSynchronizedPrintf(comm,
        "[feti user nullspace test] rank: %d, A_size: %d x %d, R_size: %d x %d, "
        "norm_A: %.16e, norm_R: %.16e, norm_A_R: %.16e, rel_A_R: %.16e, "
        "rank_estimate: %d, rank_tol: %.16e, gram_eigenvalues:",
        rank, (int)Am, (int)An, (int)Rm, (int)Rn,
        (double)norm_A, (double)norm_R, (double)norm_AR, (double)residual_rel,
        (int)rank_estimate, (double)rank_tol));
    for (PetscInt i = 0; i < Rn; ++i)
        PetscCall(PetscSynchronizedPrintf(comm, " %.16e", (double)eigenvalues[i]));
    PetscCall(PetscSynchronizedPrintf(comm, "\n"));
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));

    PetscCall(PetscFree(eigenvalues));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode print_feti_local_hessian_graph_components(MPI_Comm comm, Mat local_hessian)
{
    PetscBool enabled = PETSC_FALSE;
    PetscReal drop_tol = 0.0;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_local_hessian_graph",
                                  &enabled, NULL));
    if (!enabled) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(PetscOptionsGetReal(NULL, NULL,
                                  "-flow123d_feti_local_hessian_graph_drop_tol",
                                  &drop_tol, NULL));

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    PetscInt nrows, ncols;
    PetscCall(MatGetSize(local_hessian, &nrows, &ncols));
    PetscCheck(nrows == ncols, comm, PETSC_ERR_ARG_SIZ,
               "Local FETI Hessian graph diagnostic needs a square matrix.");

    std::vector<PetscInt> parent(nrows);
    std::vector<PetscInt> comp_size(nrows, 1);
    for (PetscInt i = 0; i < nrows; ++i) parent[i] = i;

    auto find_root = [&parent](PetscInt i) {
        PetscInt root = i;
        while (parent[root] != root) root = parent[root];
        while (parent[i] != i) {
            PetscInt next = parent[i];
            parent[i] = root;
            i = next;
        }
        return root;
    };

    auto unite = [&parent, &comp_size, &find_root](PetscInt a, PetscInt b) {
        PetscInt ra = find_root(a);
        PetscInt rb = find_root(b);
        if (ra == rb) return;
        if (comp_size[ra] < comp_size[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        comp_size[ra] += comp_size[rb];
    };

    for (PetscInt i = 0; i < nrows; ++i) {
        PetscInt ncols_row = 0;
        const PetscInt *cols = NULL;
        const PetscScalar *vals = NULL;
        PetscCall(MatGetRow(local_hessian, i, &ncols_row, &cols, &vals));
        for (PetscInt j = 0; j < ncols_row; ++j) {
            if (cols[j] == i) continue;
            if (cols[j] < 0 || cols[j] >= nrows) continue;
            if (PetscAbsScalar(vals[j]) <= drop_tol) continue;
            unite(i, cols[j]);
        }
        PetscCall(MatRestoreRow(local_hessian, i, &ncols_row, &cols, &vals));
    }

    std::unordered_map<PetscInt, PetscInt> sizes;
    for (PetscInt i = 0; i < nrows; ++i)
        ++sizes[find_root(i)];

    std::vector<PetscInt> sorted_sizes;
    sorted_sizes.reserve(sizes.size());
    for (const auto &item : sizes) sorted_sizes.push_back(item.second);
    std::sort(sorted_sizes.begin(), sorted_sizes.end(), std::greater<PetscInt>());

    PetscCall(PetscSynchronizedPrintf(comm,
        "[feti local Hessian graph] rank: %d, size: %d, components: %d, "
        "drop_tol: %.16e, component_sizes:",
        rank, (int)nrows, (int)sorted_sizes.size(), (double)drop_tol));
    for (PetscInt size : sorted_sizes)
        PetscCall(PetscSynchronizedPrintf(comm, " %d", (int)size));
    PetscCall(PetscSynchronizedPrintf(comm, "\n"));
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode print_feti_local_hessian_spectrum(MPI_Comm comm, Mat local_hessian,
                                                 IS reduced_rows, PetscInt owned_local_size)
{
    PetscBool enabled = PETSC_FALSE;
    PetscReal zero_rel_tol = 1.0e-8;
    PetscReal negative_rel_tol = 1.0e-10;
    PetscInt max_print = 16;
    PetscInt max_vector_print = 0;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_local_hessian_spectrum",
                                  &enabled, NULL));
    if (!enabled) PetscFunctionReturn(PETSC_SUCCESS);

    PetscCall(PetscOptionsGetReal(NULL, NULL,
                                  "-flow123d_feti_local_hessian_spectrum_zero_rel_tol",
                                  &zero_rel_tol, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL,
                                  "-flow123d_feti_local_hessian_spectrum_negative_rel_tol",
                                  &negative_rel_tol, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL,
                                 "-flow123d_feti_local_hessian_spectrum_max_print",
                                 &max_print, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL,
                                 "-flow123d_feti_local_hessian_spectrum_vector_split_max_print",
                                 &max_vector_print, NULL));

    PetscInt m, n;
    PetscCall(MatGetSize(local_hessian, &m, &n));
    PetscCheck(m == n, comm, PETSC_ERR_ARG_SIZ,
               "Local FETI Hessian is not square.");

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    if (n == 0) {
        PetscCall(PetscSynchronizedPrintf(comm,
            "[feti local Hessian spectrum] rank: %d, size: 0\n", rank));
        PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    Mat dense = NULL;
    PetscCall(MatConvert(local_hessian, MATSEQDENSE, MAT_INITIAL_MATRIX, &dense));

    PetscScalar *array = NULL;
    PetscCall(MatDenseGetArray(dense, &array));

    PetscBLASInt bn, lda, info;
    PetscCall(PetscBLASIntCast(n, &bn));
    lda = bn;

    PetscReal *eigenvalues = NULL;
    PetscCall(PetscMalloc1(n, &eigenvalues));

    PetscBLASInt lwork = -1;
    PetscScalar work_query;
    const char jobz = max_vector_print > 0 ? 'V' : 'N';
    LAPACKsyev_(&jobz, "U", &bn, array, &lda, eigenvalues,
                &work_query, &lwork, &info);
    PetscCheck(info == 0, comm, PETSC_ERR_LIB,
               "LAPACKsyev workspace query failed with info=%d", (int)info);

    lwork = static_cast<PetscBLASInt>(PetscRealPart(work_query));
    PetscScalar *work = NULL;
    PetscCall(PetscMalloc1(lwork, &work));

    LAPACKsyev_(&jobz, "U", &bn, array, &lda, eigenvalues,
                work, &lwork, &info);

    PetscCall(PetscFree(work));

    PetscCheck(info == 0, comm, PETSC_ERR_LIB,
               "LAPACKsyev failed with info=%d", (int)info);

    PetscReal max_abs = 0.0;
    for (PetscInt i = 0; i < n; ++i)
        max_abs = std::max(max_abs, PetscAbsReal(eigenvalues[i]));

    const PetscReal scale = max_abs > 0.0 ? max_abs : 1.0;
    const PetscReal zero_tol = zero_rel_tol * scale;
    const PetscReal negative_tol = negative_rel_tol * scale;

    PetscInt near_zero = 0;
    PetscInt negative = 0;
    for (PetscInt i = 0; i < n; ++i) {
        if (PetscAbsReal(eigenvalues[i]) <= zero_tol) ++near_zero;
        if (eigenvalues[i] < -negative_tol) ++negative;
    }

    PetscCall(PetscSynchronizedPrintf(comm,
        "[feti local Hessian spectrum] rank: %d, size: %d, lambda_min: %.16e, "
        "lambda_max: %.16e, max_abs: %.16e, near_zero: %d, negative: %d, "
        "zero_tol: %.16e, negative_tol: %.16e\n",
        rank, (int)n, (double)eigenvalues[0], (double)eigenvalues[n - 1],
        (double)max_abs, (int)near_zero, (int)negative,
        (double)zero_tol, (double)negative_tol));

    const PetscInt nprint = std::max<PetscInt>(0, std::min(n, max_print));
    PetscCall(PetscSynchronizedPrintf(comm,
        "[feti local Hessian spectrum] rank: %d, first_eigenvalues:", rank));
    for (PetscInt i = 0; i < nprint; ++i)
        PetscCall(PetscSynchronizedPrintf(comm, " %.16e", (double)eigenvalues[i]));
    if (nprint < n)
        PetscCall(PetscSynchronizedPrintf(comm, " ..."));
    PetscCall(PetscSynchronizedPrintf(comm, "\n"));

    if (max_vector_print > 0 && reduced_rows != NULL) {
        const PetscInt *rows = NULL;
        PetscInt nrows_reduced = 0;
        PetscCall(ISGetLocalSize(reduced_rows, &nrows_reduced));
        PetscCheck(nrows_reduced == n, comm, PETSC_ERR_ARG_SIZ,
                   "Reduced row IS size is not compatible with local Hessian size.");
        PetscCall(ISGetIndices(reduced_rows, &rows));

        PetscInt printed = 0;
        for (PetscInt ev = 0; ev < n && printed < max_vector_print; ++ev) {
            if (PetscAbsReal(eigenvalues[ev]) > zero_tol) continue;

            PetscReal owned_norm2 = 0.0;
            PetscReal ghost_norm2 = 0.0;
            PetscReal max_abs = 0.0;
            PetscInt max_row = -1;
            for (PetscInt row = 0; row < n; ++row) {
                const PetscScalar value = array[ev*n + row];
                const PetscReal abs_value = PetscAbsScalar(value);
                const PetscReal value_norm2 = abs_value * abs_value;
                if (rows[row] < owned_local_size) owned_norm2 += value_norm2;
                else ghost_norm2 += value_norm2;
                if (abs_value > max_abs) {
                    max_abs = abs_value;
                    max_row = rows[row];
                }
            }

            PetscCall(PetscSynchronizedPrintf(comm,
                "[feti local Hessian near-zero vector] rank: %d, eigen_index: %d, "
                "eigenvalue: %.16e, owned_norm: %.16e, ghost_norm: %.16e, "
                "max_abs: %.16e, max_original_local_row: %d\n",
                rank, (int)ev, (double)eigenvalues[ev],
                (double)std::sqrt(owned_norm2), (double)std::sqrt(ghost_norm2),
                (double)max_abs, (int)max_row));
            ++printed;
        }
        PetscCall(ISRestoreIndices(reduced_rows, &rows));
    }
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));

    PetscCall(MatDenseRestoreArray(dense, &array));
    PetscCall(MatDestroy(&dense));
    PetscCall(PetscFree(eigenvalues));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode print_feti_permon_nullspace_dimension(MPI_Comm comm, QP qp)
{
    PetscBool enabled = PETSC_FALSE;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_permon_nullspace_dimension",
                                  &enabled, NULL));
    if (!enabled) PetscFunctionReturn(PETSC_SUCCESS);

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    PetscInt chain_index = 0;
    QP current = qp;
    while (current) {
        Mat R = NULL;
        PetscErrorCode (*transform)(QP) = NULL;
        PetscCall(QPGetOperatorNullSpace(current, &R));
        PetscCall(QPGetTransform(current, &transform));

        const char *transform_name = "root";
        if (transform == (PetscErrorCode (*)(QP))QPTScale) transform_name = "QPTScale";
        else if (transform == (PetscErrorCode (*)(QP))QPTEnforceEqByProjector) transform_name = "QPTEnforceEqByProjector";
        else if (transform == (PetscErrorCode (*)(QP))QPTHomogenizeEq) transform_name = "QPTHomogenizeEq";
        else if (transform) transform_name = "other";

        if (!R) {
            PetscCall(PetscSynchronizedPrintf(comm,
                "[feti PERMON nullspace] rank: %d, qp_chain_index: %d, "
                "transform: %s, R: null, dim: 0\n",
                rank, (int)chain_index, transform_name));
        } else {
            PetscInt m_local, n_local, m_global, n_global;
            PetscCall(MatGetLocalSize(R, &m_local, &n_local));
            PetscCall(MatGetSize(R, &m_global, &n_global));

            PetscCall(PetscSynchronizedPrintf(comm,
                "[feti PERMON nullspace] rank: %d, qp_chain_index: %d, "
                "transform: %s, R_local_size: %d x %d, "
                "R_global_size: %d x %d, dim: %d\n",
                rank, (int)chain_index, transform_name,
                (int)m_local, (int)n_local, (int)m_global, (int)n_global,
                (int)n_global));
        }

        PetscCall(QPGetChild(current, &current));
        ++chain_index;
    }
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscScalar dense_row_dot(const std::vector<PetscScalar> &left,
                          const std::vector<PetscScalar> &right)
{
    PetscScalar dot = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
        dot += PetscConj(left[i]) * right[i];
    return dot;
}

PetscReal dense_row_norm2(const std::vector<PetscScalar> &row)
{
    PetscReal norm2 = 0.0;
    for (PetscScalar value : row) {
        const PetscReal abs_value = PetscAbsScalar(value);
        norm2 += abs_value * abs_value;
    }
    return norm2;
}

void dense_row_subtract_scaled(std::vector<PetscScalar> &row,
                               PetscScalar scale,
                               const std::vector<PetscScalar> &basis_row)
{
    for (std::size_t i = 0; i < row.size(); ++i)
        row[i] -= scale * basis_row[i];
}

PetscErrorCode gather_dense_matrix_rows(Mat matrix,
                                        std::vector<PetscScalar> &dense_rows,
                                        PetscInt *global_rows,
                                        PetscInt *global_cols)
{
    Mat matrix_aij = NULL;
    PetscInt m_global, n_global, rstart, rend;
    MPI_Comm comm;
    PetscMPIInt comm_size;

    PetscFunctionBegin;
    PetscCall(PetscObjectGetComm((PetscObject)matrix, &comm));
    PetscCallMPI(MPI_Comm_size(comm, &comm_size));

    PetscCall(MatConvert(matrix, MATAIJ, MAT_INITIAL_MATRIX, &matrix_aij));
    PetscCall(MatGetSize(matrix_aij, &m_global, &n_global));
    PetscCall(MatGetOwnershipRange(matrix_aij, &rstart, &rend));

    PetscCheck(n_global <= static_cast<PetscInt>(std::numeric_limits<PetscMPIInt>::max()),
               comm, PETSC_ERR_SUP,
               "FETI equality-rank filter matrix is too wide for MPI gather diagnostics.");

    const PetscInt local_rows = rend - rstart;
    PetscCheck(local_rows == 0 ||
               n_global <= std::numeric_limits<PetscInt>::max() / local_rows,
               comm, PETSC_ERR_SUP,
               "FETI equality-rank filter matrix block is too large for dense row gathering.");
    const PetscInt local_entries = local_rows * n_global;
    PetscMPIInt local_count;
    PetscCall(PetscMPIIntCast(local_entries, &local_count));

    std::vector<PetscScalar> local_dense(local_entries, 0.0);
    for (PetscInt row = rstart; row < rend; ++row) {
        PetscInt ncols = 0;
        const PetscInt *cols = NULL;
        const PetscScalar *vals = NULL;

        PetscCall(MatGetRow(matrix_aij, row, &ncols, &cols, &vals));
        for (PetscInt j = 0; j < ncols; ++j)
            local_dense[(row - rstart) * n_global + cols[j]] = vals[j];
        PetscCall(MatRestoreRow(matrix_aij, row, &ncols, &cols, &vals));
    }

    std::vector<PetscMPIInt> recv_counts(comm_size, 0);
    std::vector<PetscMPIInt> recv_displs(comm_size, 0);
    PetscCallMPI(MPI_Allgather(&local_count, 1, MPI_INT,
                               recv_counts.data(), 1, MPI_INT, comm));

    PetscMPIInt total_count = 0;
    for (PetscMPIInt rank = 0; rank < comm_size; ++rank) {
        recv_displs[rank] = total_count;
        PetscCheck(total_count <= std::numeric_limits<PetscMPIInt>::max() - recv_counts[rank],
                   comm, PETSC_ERR_SUP,
                   "FETI equality-rank filter matrix is too large for MPI dense row gathering.");
        total_count += recv_counts[rank];
    }

    PetscCheck(m_global == 0 ||
               static_cast<std::size_t>(n_global) <=
               std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(m_global),
               comm, PETSC_ERR_SUP,
               "FETI equality-rank filter matrix is too large for dense row storage.");
    dense_rows.assign(static_cast<std::size_t>(m_global) * n_global, 0.0);
    PetscCallMPI(MPI_Allgatherv(local_dense.data(), local_count, MPIU_SCALAR,
                                dense_rows.data(), recv_counts.data(), recv_displs.data(),
                                MPIU_SCALAR, comm));

    PetscCall(MatDestroy(&matrix_aij));

    if (global_rows) *global_rows = m_global;
    if (global_cols) *global_cols = n_global;
    PetscFunctionReturn(PETSC_SUCCESS);
}

std::vector<PetscScalar> dense_matrix_row(const std::vector<PetscScalar> &dense_rows,
                                          PetscInt ncols,
                                          PetscInt row)
{
    auto begin = dense_rows.begin() + static_cast<std::size_t>(row) * ncols;
    return std::vector<PetscScalar>(begin, begin + ncols);
}

bool add_row_to_rank_basis(const std::vector<PetscScalar> &input_row,
                           PetscReal row_norm2,
                           PetscReal tolerance,
                           PetscReal scale,
                           std::vector<std::vector<PetscScalar>> &basis,
                           PetscReal *residual_norm2)
{
    std::vector<PetscScalar> residual = input_row;

    // Two-pass modified Gram-Schmidt is sufficient here because this optional
    // filter is limited to small equality systems and only decides which
    // appended gluing rows are numerically dependent on the preceding rows.
    for (const auto &basis_row : basis)
        dense_row_subtract_scaled(residual, dense_row_dot(basis_row, residual), basis_row);
    for (const auto &basis_row : basis)
        dense_row_subtract_scaled(residual, dense_row_dot(basis_row, residual), basis_row);

    PetscReal res2 = dense_row_norm2(residual);
    if (residual_norm2) *residual_norm2 = res2;

    const PetscReal threshold = tolerance * std::max(scale, row_norm2);
    if (res2 <= threshold) return false;

    const PetscReal inv_norm = 1.0 / std::sqrt(res2);
    for (PetscScalar &value : residual) value *= inv_norm;
    basis.push_back(std::move(residual));
    return true;
}

PetscErrorCode filter_dependent_feti_gluing_rows(MPI_Comm comm, QP qp)
{
    PetscBool enabled = PETSC_FALSE;
    PetscReal tolerance = 1e-10;
    PetscInt max_rows = 5000;
    Mat BE = NULL;
    Vec cE = NULL;
    PetscBool is_nested = PETSC_FALSE;
    PetscInt nrow_blocks = 0, ncol_blocks = 0;

    PetscFunctionBegin;
    PetscCall(PetscOptionsGetBool(NULL, NULL,
                                  "-flow123d_feti_filter_dependent_gluing",
                                  &enabled, NULL));
    if (!enabled) PetscFunctionReturn(PETSC_SUCCESS);

    PetscCall(PetscOptionsGetReal(NULL, NULL,
                                  "-flow123d_feti_filter_dependent_gluing_tol",
                                  &tolerance, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL,
                                 "-flow123d_feti_filter_dependent_gluing_max_rows",
                                 &max_rows, NULL));

    // Flow123d adds Dirichlet equalities before QPFetiSetUp(). PERMON then
    // appends the gluing matrix as the last equality block. Keep the original
    // rows and remove only gluing rows whose residual against their row span is
    // below the requested tolerance.
    PetscCall(QPGetEq(qp, &BE, &cE));
    if (!BE) PetscFunctionReturn(PETSC_SUCCESS);

    PetscCall(PetscObjectTypeCompareAny((PetscObject)BE, &is_nested,
                                        MATNEST, MATNESTPERMON, ""));
    if (!is_nested) {
        WarningOut() << "Skipping FETI gluing rank filter: equality matrix is not nested.\n";
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscCall(MatNestGetSize(BE, &nrow_blocks, &ncol_blocks));
    if (ncol_blocks != 1 || nrow_blocks < 2) {
        WarningOut() << "Skipping FETI gluing rank filter: unexpected equality matrix nest shape.\n";
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    std::vector<Mat> blocks(nrow_blocks, NULL);
    std::vector<Vec> rhs_blocks(nrow_blocks, NULL);
    PetscBool has_nested_rhs = PETSC_FALSE;
    if (cE)
        PetscCall(PetscObjectTypeCompare((PetscObject)cE, VECNEST, &has_nested_rhs));
    if (cE && !has_nested_rhs) {
        WarningOut() << "Skipping FETI gluing rank filter: equality RHS is not nested.\n";
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    for (PetscInt block = 0; block < nrow_blocks; ++block) {
        PetscCall(MatNestGetSubMat(BE, block, 0, &blocks[block]));
        PetscCheck(blocks[block], comm, PETSC_ERR_ARG_WRONG,
                   "Nested equality matrix contains a null block.");
        PetscCall(PetscObjectReference((PetscObject)blocks[block]));
        if (has_nested_rhs) {
            PetscCall(VecNestGetSubVec(cE, block, &rhs_blocks[block]));
            PetscCall(PetscObjectReference((PetscObject)rhs_blocks[block]));
        }
    }

    const PetscInt gluing_block = nrow_blocks - 1;
    PetscInt total_rows = 0;
    PetscInt ncols_reference = -1;
    for (PetscInt block = 0; block < nrow_blocks; ++block) {
        PetscInt rows, cols;
        PetscCall(MatGetSize(blocks[block], &rows, &cols));
        total_rows += rows;
        if (ncols_reference < 0) ncols_reference = cols;
        PetscCheck(cols == ncols_reference, comm, PETSC_ERR_ARG_SIZ,
                   "Nested equality matrix blocks have inconsistent column counts.");
    }

    if (total_rows > max_rows) {
        WarningOut() << "Skipping FETI gluing rank filter: equality row count "
                     << total_rows << " exceeds limit " << max_rows << ".\n";
        for (PetscInt block = 0; block < nrow_blocks; ++block) {
            if (rhs_blocks[block]) PetscCall(VecDestroy(&rhs_blocks[block]));
            PetscCall(MatDestroy(&blocks[block]));
        }
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    std::vector<std::vector<PetscScalar>> basis;
    basis.reserve(total_rows);
    PetscReal scale = 0.0;
    PetscInt fixed_rows = 0;
    PetscInt fixed_rank = 0;
    PetscInt fixed_dependent = 0;

    for (PetscInt block = 0; block < gluing_block; ++block) {
        std::vector<PetscScalar> dense_rows;
        PetscInt rows, cols;
        PetscCall(gather_dense_matrix_rows(blocks[block], dense_rows, &rows, &cols));
        for (PetscInt row = 0; row < rows; ++row) {
            std::vector<PetscScalar> row_values = dense_matrix_row(dense_rows, cols, row);
            PetscReal row_norm2 = dense_row_norm2(row_values);
            scale = std::max(scale, row_norm2);
            PetscReal residual_norm2;
            if (add_row_to_rank_basis(row_values, row_norm2, tolerance, scale, basis, &residual_norm2))
                ++fixed_rank;
            else
                ++fixed_dependent;
            ++fixed_rows;
        }
    }

    std::vector<PetscScalar> gluing_dense_rows;
    PetscInt gluing_rows, gluing_cols;
    PetscCall(gather_dense_matrix_rows(blocks[gluing_block], gluing_dense_rows,
                                       &gluing_rows, &gluing_cols));

    std::vector<PetscInt> keep_gluing_rows;
    keep_gluing_rows.reserve(gluing_rows);
    PetscReal min_accepted_residual = PETSC_MAX_REAL;
    PetscReal max_rejected_residual = 0.0;
    for (PetscInt row = 0; row < gluing_rows; ++row) {
        std::vector<PetscScalar> row_values = dense_matrix_row(gluing_dense_rows, gluing_cols, row);
        PetscReal row_norm2 = dense_row_norm2(row_values);
        scale = std::max(scale, row_norm2);
        PetscReal residual_norm2;
        if (add_row_to_rank_basis(row_values, row_norm2, tolerance, scale, basis, &residual_norm2)) {
            keep_gluing_rows.push_back(row);
            min_accepted_residual = std::min(min_accepted_residual, residual_norm2);
        } else {
            max_rejected_residual = std::max(max_rejected_residual, residual_norm2);
        }
    }

    PetscInt rstart, rend, cstart, cend;
    PetscCall(MatGetOwnershipRange(blocks[gluing_block], &rstart, &rend));
    PetscCall(MatGetOwnershipRangeColumn(blocks[gluing_block], &cstart, &cend));

    PetscCheck(keep_gluing_rows.size() <=
               static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()),
               comm, PETSC_ERR_SUP,
               "Too many globally kept FETI gluing rows.");
    std::vector<PetscInt> local_keep_rows;
    std::vector<PetscInt> local_keep_new_rows;
    for (PetscInt new_row = 0; new_row < static_cast<PetscInt>(keep_gluing_rows.size()); ++new_row) {
        const PetscInt source_row = keep_gluing_rows[new_row];
        if (source_row >= rstart && source_row < rend) {
            local_keep_rows.push_back(source_row);
            local_keep_new_rows.push_back(new_row);
        }
    }

    Mat filtered_gluing = NULL;
    PetscCheck(local_keep_rows.size() <=
               static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()),
               comm, PETSC_ERR_SUP,
               "Too many locally kept FETI gluing rows.");
    const PetscInt n_local_keep_rows = static_cast<PetscInt>(local_keep_rows.size());
    const PetscInt n_global_keep_rows = static_cast<PetscInt>(keep_gluing_rows.size());
    std::vector<PetscInt> diagonal_nnz(n_local_keep_rows, 0);
    std::vector<PetscInt> offdiagonal_nnz(n_local_keep_rows, 0);
    for (PetscInt local_row = 0; local_row < n_local_keep_rows; ++local_row) {
        const PetscInt source_row = local_keep_rows[local_row];
        for (PetscInt col = 0; col < gluing_cols; ++col) {
            const PetscScalar value =
                gluing_dense_rows[static_cast<std::size_t>(source_row) * gluing_cols + col];
            if (PetscAbsScalar(value) == 0.0) continue;
            if (col >= cstart && col < cend)
                ++diagonal_nnz[local_row];
            else
                ++offdiagonal_nnz[local_row];
        }
    }

    PetscCall(MatCreateAIJ(comm, n_local_keep_rows, cend - cstart,
                           n_global_keep_rows, gluing_cols,
                           0, diagonal_nnz.data(), 0, offdiagonal_nnz.data(),
                           &filtered_gluing));
    PetscInt filtered_rstart, filtered_rend;
    PetscCall(MatGetOwnershipRange(filtered_gluing, &filtered_rstart, &filtered_rend));
    if (n_local_keep_rows > 0) {
        PetscCheck(filtered_rstart == local_keep_new_rows.front() &&
                   filtered_rend == local_keep_new_rows.back() + 1,
                   comm, PETSC_ERR_ARG_WRONGSTATE,
                   "Filtered FETI gluing row ownership is inconsistent.");
    }
    for (PetscInt local_row = 0; local_row < n_local_keep_rows; ++local_row) {
        const PetscInt source_row = local_keep_rows[local_row];
        const PetscInt new_row = local_keep_new_rows[local_row];
        std::vector<PetscInt> cols;
        std::vector<PetscScalar> vals;
        cols.reserve(diagonal_nnz[local_row] + offdiagonal_nnz[local_row]);
        vals.reserve(diagonal_nnz[local_row] + offdiagonal_nnz[local_row]);
        for (PetscInt col = 0; col < gluing_cols; ++col) {
            const PetscScalar value =
                gluing_dense_rows[static_cast<std::size_t>(source_row) * gluing_cols + col];
            if (PetscAbsScalar(value) == 0.0) continue;
            cols.push_back(col);
            vals.push_back(value);
        }
        PetscCall(MatSetValues(filtered_gluing, 1, &new_row,
                               static_cast<PetscInt>(cols.size()),
                               cols.empty() ? NULL : cols.data(),
                               vals.empty() ? NULL : vals.data(),
                               INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(filtered_gluing, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(filtered_gluing, MAT_FINAL_ASSEMBLY));
    PetscCall(PetscObjectSetName((PetscObject)filtered_gluing, "Bg_filtered"));

    PetscCall(QPSetEq(qp, NULL, NULL));
    PetscCall(QPSetEq(qp, blocks[0], rhs_blocks[0]));
    for (PetscInt block = 1; block < gluing_block; ++block)
        PetscCall(QPAddEq(qp, blocks[block], rhs_blocks[block]));
    PetscCall(QPAddEq(qp, filtered_gluing, NULL));

    if (min_accepted_residual == PETSC_MAX_REAL) min_accepted_residual = 0.0;
    PetscCall(PetscPrintf(comm,
        "[flow123d feti gluing rank filter] fixed_rows: %d, fixed_rank: %d, "
        "fixed_dependent: %d, gluing_rows_before: %d, gluing_rows_kept: %d, "
        "gluing_rows_removed: %d, rank_estimate_after: %d, tolerance: %.16e, "
        "min_accepted_residual2: %.16e, max_rejected_residual2: %.16e\n",
        (int)fixed_rows, (int)fixed_rank, (int)fixed_dependent,
        (int)gluing_rows, (int)keep_gluing_rows.size(),
        (int)(gluing_rows - static_cast<PetscInt>(keep_gluing_rows.size())),
        (int)basis.size(), (double)tolerance,
        (double)min_accepted_residual, (double)max_rejected_residual));

    PetscCall(MatDestroy(&filtered_gluing));
    for (PetscInt block = 0; block < nrow_blocks; ++block) {
        if (rhs_blocks[block]) PetscCall(VecDestroy(&rhs_blocks[block]));
        PetscCall(MatDestroy(&blocks[block]));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode compute_qp_lagrangian_residual_norm(QP qp, PetscReal *norm)
{
    Vec solution = NULL;
    Vec residual = NULL;
    char *residual_name = NULL;
    PetscBool invalid_residual = PETSC_FALSE;

    PetscFunctionBegin;
    PetscCheck(qp, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null QP for residual computation.");
    PetscCheck(norm, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null residual norm pointer.");

    *norm = -1.0;
    PetscCall(QPComputeMissingEqMultiplier(qp));
    PetscCall(QPComputeMissingBoxMultipliers(qp));
    PetscCall(QPGetSolutionVector(qp, &solution));
    PetscCall(VecDuplicate(solution, &residual));
    PetscCall(QPComputeLagrangianGradient(qp, solution, residual, &residual_name));
    PetscCall(VecIsInvalidated(residual, &invalid_residual));
    if (!invalid_residual)
        PetscCall(VecNorm(residual, NORM_2, norm));
    PetscCall(PetscFree(residual_name));
    PetscCall(VecDestroy(&residual));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

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
          matrix_eq_(NULL),
          eq_(NULL),
          operator_nullspace_(NULL),
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
	: LinSys(other), params_(other.params_), l2g_(other.l2g_), solution_precision_(other.solution_precision_),
      operator_nullspace_(other.operator_nullspace_)
{
    MatCopy(other.matrix_, matrix_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.rhs_, rhs_);
	VecCopy(other.on_vec_, on_vec_);
	VecCopy(other.off_vec_, off_vec_);

	MatCopy(other.matrix_ineq_, matrix_ineq_, DIFFERENT_NONZERO_PATTERN);
	VecCopy(other.ineq_, ineq_);
    MatCopy(other.matrix_eq_, matrix_eq_, DIFFERENT_NONZERO_PATTERN);
    VecCopy(other.eq_, eq_);
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


void LinSys_PERMON::set_equality(Mat matrix_eq, Vec eq)
{
    matrix_eq_ = matrix_eq;
    eq_ = eq;
}

void LinSys_PERMON::set_operator_nullspace(Mat operator_nullspace)
{
    // TODO ref count?
    operator_nullspace_ = operator_nullspace;
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
    // const char *petsc_dflt_opt;
    int nits;
    
    // -mat_no_inode ... inodes are usefull only for
    //  vector problems e.g. MH without Schur complement reduction
    
    /* Comment to PETSc options:
     * 
     * -ksp_diagonal_scale scales the matrix before solution, while -ksp_diagonal_scale_fix just fixes the scaling after solution
     * -pc_asm_type basic enforces classical Schwartz method, which seems more stable for positive definite systems.
     *                    The default 'restricted' probably violates s.p.d. structure, many tests fail.
     */
    // if (rows_ds_->np() > 1) {
    //     // parallel setting
    //    if (this->is_positive_definite())
    //        petsc_dflt_opt="-ksp_type cg -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_type asm -pc_asm_type basic -pc_asm_overlap 4 -sub_pc_type icc -sub_pc_factor_levels 3  -sub_pc_factor_fill 6.0";
    //        //petsc_dflt_opt="-ksp_type bcgs -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3  -sub_pc_factor_fill 6.0";
    //    else
    //        petsc_dflt_opt="-ksp_type bcgs -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_type asm -pc_asm_overlap 4 -sub_pc_type ilu -sub_pc_factor_levels 3 -sub_pc_factor_fill 6.0";
    
    // }
    // else {
    //     // serial setting
    //    if (this->is_positive_definite())
    //        petsc_dflt_opt="-ksp_type cg -pc_type icc  -pc_factor_levels 3 -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    // 	   //petsc_dflt_opt="-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    //    else
    //        petsc_dflt_opt="-ksp_type bcgs -pc_type ilu -pc_factor_levels 5 -ksp_diagonal_scale -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
    // }

    // if (params_ == "") params_ = petsc_dflt_opt;
    LogOut().fmt("inserting petsc options: {}\n",params_.c_str());
    
    // now takes an optional PetscOptions object as the first argument
    // value NULL will preserve previous behaviour previous behavior.
    PetscOptionsInsertString(NULL, params_.c_str()); // overwrites previous options values
    
    MatSetOption( matrix_, MAT_USE_INODES, PETSC_FALSE );

    if (use_feti_ && (rows_ds_->np() == 1)) {
        WarningOut() << "FETI can be only used on multiple processes - switching to standard QP solver.";
        use_feti_ = false;
    }

    QP primal_qp = NULL;

    chkerr(QPCreate(comm_, &system));
    chkerr(QPSetOptionsPrefix(system,"permon_")); // avoid clash on PC objects from hydro PETSc solver

    chkerr(QPSetRhs(system, rhs_));
    chkerr(QPSetInitialVector(system, solution_));
    {
        Vec qp_solution = NULL;
        PetscReal qp_solution_norm;
        chkerr(QPGetSolutionVector(system, &qp_solution));
        chkerr(VecNorm(qp_solution, NORM_2, &qp_solution_norm));
        if (!std::isfinite(qp_solution_norm)) {
            WarningOut() << "PERMON QP initial vector contains NaN or Inf values. "
                         << "Resetting the initial vector to zero.\n";
            chkerr(VecSet(qp_solution, 0.0));
        }
    }

    if (use_feti_) {

        // create new matrix without zero rows
        Mat Afixed;
        Mat Aloc, Aloc_feti;
        IS isnz, isnz_tmp = NULL;
        ISLocalToGlobalMapping l2g, mapping;
        const PetscInt *l2g_arr, *nz_arr;
        PetscInt *l2g_feti_arr;
        PetscInt nmap, nrows, k;
        
        chkerr(MatISGetLocalMat(matrix_, &Aloc));
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
        chkerr(print_feti_local_hessian_graph_components(comm_, Aloc_feti));
        chkerr(print_feti_local_hessian_spectrum(comm_, Aloc_feti, isnz, rows_ds_->lsize()));
        chkerr(MatISRestoreLocalMat(matrix_, &Aloc));
        chkerr(ISGetIndices(isnz, &nz_arr));
        chkerr(ISGetLocalSize(isnz, &nrows));
        chkerr(MatISGetLocalToGlobalMapping(matrix_, &l2g, NULL));
        chkerr(ISLocalToGlobalMappingGetIndices(l2g, &l2g_arr));
        chkerr(ISLocalToGlobalMappingGetSize(l2g, &nmap));
        chkerr(PetscMalloc1(nrows, &l2g_feti_arr));
        for (k=0; k<nrows; k++) {
            ASSERT_LT(static_cast<unsigned int>(nz_arr[k]), static_cast<unsigned int>(nmap))
                .error("Reduced FETI local row index out of bounds of local-to-global mapping.\n");
            l2g_feti_arr[k] = l2g_arr[nz_arr[k]];
        }
        chkerr(ISRestoreIndices(isnz, &nz_arr));
        chkerr(ISLocalToGlobalMappingCreate(comm_, 1, nrows, l2g_feti_arr, PETSC_OWN_POINTER, &mapping));
        chkerr(ISLocalToGlobalMappingRestoreIndices(l2g, &l2g_arr));
        
        chkerr(MatCreateIS(PETSC_COMM_WORLD, 1, rows_ds_->lsize(), rows_ds_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                            mapping, mapping, &Afixed));
        chkerr(MatISSetLocalMat(Afixed, Aloc_feti));
        chkerr(MatAssemblyBegin(Afixed, MAT_FINAL_ASSEMBLY));
        chkerr(MatAssemblyEnd(Afixed, MAT_FINAL_ASSEMBLY));

        // set QP matrix
        chkerr(QPSetOperator(system, Afixed));
        chkerr(MatDestroy(&Afixed));

        chkerr(QPTMatISToBlockDiag(system));
        chkerr(QPGetChild(system, &system)); // now system is decomposed

        Mat Adecomposed = NULL;
        {
            QPGetOperator(system, &Adecomposed);
            if (matrix_ineq_) {
                Mat matrix_ineq_decomposed = NULL;
                chkerr(convert_mat_is_to_feti_decomposed(matrix_ineq_, isnz, Adecomposed, &matrix_ineq_decomposed));
                chkerr(QPSetIneq(system, matrix_ineq_decomposed, ineq_));
                chkerr(MatDestroy(&matrix_ineq_decomposed));
                chkerr(PetscOptionsInsertString(NULL, "-qpt_dualize_B_nest_extension"));
            }
            if (matrix_eq_) {
                Mat matrix_eq_decomposed = NULL;
                chkerr(convert_mat_is_to_feti_decomposed(matrix_eq_, isnz, Adecomposed, &matrix_eq_decomposed));
                chkerr(QPSetEq(system, matrix_eq_decomposed, eq_));
                chkerr(MatDestroy(&matrix_eq_decomposed));
                chkerr(PetscOptionsInsertString(NULL, "-qpt_dualize_B_nest_extension"));
            }
        }
        if (operator_nullspace_) {
            Mat Rloc_reduced = NULL;
            Mat nullspace_decomposed = NULL;
            chkerr(create_reduced_nullspace_local(operator_nullspace_, isnz, &Rloc_reduced));
            chkerr(print_feti_user_nullspace_test(comm_, Aloc_feti, Rloc_reduced));
            chkerr(convert_nullspace_to_feti_blockdiag(Rloc_reduced, &nullspace_decomposed));
            chkerr(print_feti_user_nullspace_dimension(comm_, nullspace_decomposed));
            chkerr(QPSetOperatorNullSpace(system, nullspace_decomposed));
            chkerr(MatDestroy(&nullspace_decomposed));
            chkerr(MatDestroy(&Rloc_reduced));
            chkerr(PetscOptionsInsertString(NULL, "-qpt_dualize_Kplus_mp -regularize 0"));
        }
        chkerr(MatDestroy(&Aloc_feti));
        chkerr(ISLocalToGlobalMappingDestroy(&mapping));
        chkerr(ISDestroy(&isnz));
        chkerr(QPFetiSetUp(system));
        chkerr(filter_dependent_feti_gluing_rows(comm_, system));
        chkerr(PetscOptionsInsertString(NULL, "-feti"));

        // Set/Unset additional transformations, e.g -project 0 for projector avoiding FETI
        chkerr(QPTFromOptions(system));
        chkerr(print_feti_permon_nullspace_dimension(comm_, system));
        primal_qp = system;
        chkerr(QPGetParent(system, &system));
    } else {
        // convert to MATAIJ
        Mat matrix_aij;
        chkerr(MatConvert(matrix_, MATAIJ, MAT_INITIAL_MATRIX, &matrix_aij));
        if (eq_) {
            Mat matrix_eq_aij = NULL;
            chkerr(convert_mat_is_to_aij(matrix_eq_, &matrix_eq_aij));
            chkerr(QPSetEq(system, matrix_eq_aij, eq_));
        }
        else {
            chkerr(MatSetOption(matrix_aij,MAT_SPD,PETSC_TRUE)); // avoid null space computation
        // chkerr(MatSetOption(matrix_aij,MAT_SPD_ETERNAL,PETSC_TRUE)); // possible with PETSc >= 3.18.0
        }
        if (ineq_) {
            Mat matrix_ineq_aij = NULL;
            chkerr(convert_mat_is_to_aij(matrix_ineq_, &matrix_ineq_aij));
            chkerr(QPSetIneq(system, matrix_ineq_aij, ineq_));
            chkerr(MatDestroy(&matrix_ineq_aij));
        }
        chkerr(QPSetOperator(system, matrix_aij));
        chkerr(MatDestroy(&matrix_aij));

        if (ineq_ || eq_) // dualization without FETI
            chkerr(QPTDualize(system, MAT_INV_MONOLITHIC, MAT_REG_NONE));
        if (ineq_ || eq_)
            primal_qp = system;
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
    if (primal_qp) {
        // Compute the Lagrangian residual |A*x - b + B'*lambda| on the
        // primal QP in both FETI and constrained non-FETI paths.
        chkerr(compute_qp_lagrangian_residual_norm(primal_qp, &solution_precision_));
    } else {
        MatMult(matrix_, solution_, residual_);
        VecAXPY(residual_,-1.0, rhs_);
        VecNorm(residual_, NORM_2, &solution_precision_);
    }

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
    FilePath mat_eqFileName(text + "_flow123d_matrix_eq.m",FilePath::FileType::output_file);
    FilePath eqFileName(text + "_flow123d_eq.m",FilePath::FileType::output_file);

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

    if ( matrix_eq_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)mat_eqFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        MatView( matrix_eq_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the equality matrix of LinSys is not set.\n";
    if ( eq_ != NULL ) {
        PetscViewerASCIIOpen(comm_,((string)eqFileName).c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( eq_, myViewer );
        PetscViewerDestroy(&myViewer);
    }
    else
        WarningOut() << "PetscViewer: the equality vector of LinSys is not set.\n";
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
    return solution_precision_;
}
