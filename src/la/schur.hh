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
 *
 * @file
 * @brief  Assembly explicit Schur complement for the given linear system.
 * Provides method for resolution of the full original vector of unknowns.
 *
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
 *  TODO:
 *  - make derived classes for MPI and MATIS matrices
 *  - make LinSys have its own solve method
 *  - make consistent work cycle:
 *    1) Constructor only takes original system and index set of reduced block
 *    2) solve - takes inv_a as parameter , optional parameter different_structure
 *       - when called first call form_schur without reuse of matrices, latter calls use reuse
 *       - different_structure is set -> either inv_a or orig_linsys structure changed
 *         has to call destruction of auxiliar matricies and call initial phase of form_schur
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
    SchurComplement(LinSys *orig, IS ia, Distribution *ds);
    SchurComplement(IS ia, Distribution *ds);

    /**
     * Copy constructor.
     */
    SchurComplement(SchurComplement &other);

    LinSys *get_orig_system() const {return (Orig);}
    Distribution *get_distribution() const {return (ds_);}
    Mat get_a_inv() const {return (IA);}

    void scale(double factor);
    ~SchurComplement();

    // TODO: should be at least protected
    void form_schur();
    /** Compute only right hand side.
     *  This is useful when you change only rhs of the original system.
     *  TODO: We should ask original system if the matrix has changed (using LazyDependency) and
     *  possibly call only form_rhs, then this can be protected
     */
    void form_rhs();
    void resolve();
    /// set complement object
    void set_complement(LinSys_PETSC *ls);
    /// get distribution of complement object if complement is defined
    Distribution *make_complement_distribution();
    /// create IA matrix and compute complement
    void create_inversion_matrix();

private:
    Mat IA;                     // Inverse of block A

    Mat B, Bt;                  // B and B' block (could be different from real B transpose)
    Mat xA;                     // Bt*IA*B
    Mat IAB;                    // reconstruction matrix IA * B
    int loc_size_A, loc_size_B; // loc size of the A and B block
    IS IsA, IsB;                // parallel index sets of the A and B block
    IS fullIsA,fullIsB;         // whole IsA  and IsB on each proc
    Vec RHS1, RHS2;             // A and B - part of the RHS
    Vec Sol1, Sol2;             // A and B part of solution

    SchurState state;           // object internal state
    int orig_lsize;             ///< Size of local vector part of original system

                                //                A  B     Sol1      RHS1
    LinSys *Orig;               // Original Linear System:  B' C  *  Sol2  =   RHS2
    LinSys_PETSC *Compl;        // Schur complement system: (C - B' IA B) * Sol2 = (B' * IA * RHS1 - RHS2)

    Distribution *ds_;          // Distribution of B block
} SchurComplement;

#endif /* LA_SCHUR_HH_ */
