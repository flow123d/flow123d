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
 * @file    schur.hh
 * @brief   Assembly explicit Schur complement for the given linear system.
 *          Provides method for resolution of the full original vector of unknowns.
 */

#ifndef LA_SCHUR_HH_
#define LA_SCHUR_HH_

#include <la/distribution.hh>
#include "la/linsys_PETSC.hh"
#include <petscmat.h>

struct Solver;
class LinSys;

/**
 * @brief Schur complement class for a PETSC based linear system
 *
 * Linear system consists of Matrix, inversion matrix, RHS, solution and pointer to Complement.
 * It provides methods for:
 * - set complement system
 * - create distribution of complement system
 * - create inversion matrix of system
 * - form Schur complement and rhs
 * - solve and resolve system
 *
 * Usage example:
 * @CODE
 *   SchurComplement * schur = new SchurComplement(is, &distr);
 *
 *   ... // allocation and assembly of matrix
 *
 *   LinSys_PETSC * ls = new LinSys_PETSC( schur->make_complement_distribution() );
 *   schur->set_complement( ls );
 *   schur->solve();
 * @ENDCODE
 *
 * Input record is passed to the complement system.
 */

typedef enum SchurState {
    created,    // need creation of all temporal matrixes ...
    formed,     // Schur complement ready to solve
    solved      // Resolved original Schur system
} SchurState;

typedef class SchurComplement : public LinSys_PETSC {
public:
    /**
     * Constructor
     *
     * Gets linear system with original matrix A and creates its inversion (IA matrix)
     *
     * In current implementation the index set IsA has to be continuous sequence at the beginning of the local block of indices.
     */
    SchurComplement(IS ia, Distribution *ds);

    /**
     * Copy constructor.
     */
    SchurComplement(SchurComplement &other);

    /**
     * Returns pointer to LinSys object representing the schur complement.
     */
    LinSys *get_system() const {return (Compl);}

    /**
     * Returns distribution of the original system (solved by class SchurComplement).
     */
    Distribution *get_distribution() const {return (ds_);}

    /**
     * Destructor. In particular it also delete complement linear system if it was passed in
     * through the @p set_complement() method.
     */
    ~SchurComplement();

    /** Compute only right hand side.
     *  This is useful when you change only rhs of the original system.
     *  TODO: We should ask original system if the matrix has changed (using LazyDependency) and
     *  possibly call only form_rhs, then this can be protected
     */
    void form_rhs();
    /// Set complement LinSys object.
    void set_complement(LinSys_PETSC *ls);
    /// get distribution of complement object if complement is defined
    Distribution *make_complement_distribution();
    /// get precision of solving
    double get_solution_precision() override;

    /**
     * Solve the system.
     *
     * - If matrix and/or rhs is changed the Schur complement is formed.
     * - The inner linear solver is called for the Schur complement
     * - Resolve is called to reconstruct eliminated part of the solution vector.
     */
    int solve() override;

    /**
     * Only resolve the system with current solution vector. This is necessary for nonlinear solvers.
     * - If matrix and/or rhs is changed the Schur complement is formed.
     * - Resolve is called to reconstruct eliminated part of the solution vector.
     */
    void resolve();

    /**
     * The solve or resolve must be called prior to computing the residual.
     */
    double compute_residual() override;

protected:
    /// create IA matrix
    void create_inversion_matrix();

    void form_schur();



    Mat IA;                     // Inverse of block A

    Mat B, Bt;                  // B and B' block (could be different from real B transpose)
    Mat xA;                     // Bt*IA*B
    Mat IAB;                    // reconstruction matrix IA * B
    int loc_size_A, loc_size_B; // loc size of the A and B block
    IS IsA, IsB;                // parallel index sets of the A and B block
    Vec RHS1, RHS2;             // A and B - part of the RHS
    Vec Sol1, Sol2;             // A and B part of solution

    SchurState state;           // object internal state
    int orig_lsize;             ///< Size of local vector part of original system

    LinSys_PETSC *Compl;        // Schur complement system: (C - B' IA B) * Sol2 = (B' * IA * RHS1 - RHS2)

    Distribution *ds_;          // Distribution of B block
} SchurComplement;

#endif /* LA_SCHUR_HH_ */
