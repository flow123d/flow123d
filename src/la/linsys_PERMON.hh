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
 * @file    linsys_PERMON.hh
 * @brief   PERMON QP solvers and FETI
 * @author  Jakub Kruzik
 */

#ifndef LA_LINSYS_PERMON_HH_
#define LA_LINSYS_PERMON_HH_

#include <functional>          // for unary_function
#include <string>              // for string
#include <vector>              // for vector
#include "la/linsys_PETSC.hh"  // for LinSys
#include "permonqps.h"         // for QPS and whole PERMON/PETSc stack

class Distribution;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}
namespace la {
    class BddcmlWrapper;
}

class LinSys_PERMON : public LinSys_PETSC
{

public:
	typedef LinSys FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    LinSys_PERMON(const  Distribution * rows_ds, const std::string &params = "");

    /**
     * Copy constructor.
     */
    LinSys_PERMON( LinSys_PERMON &other );

    /// Set A*x <= b using borrowed objects in original global primal indexing; caller keeps them alive.
    void set_inequality(Mat matrix_ineq, Vec ineq);

    /**
     * Add the support function max_{lower <= q <= upper} q^T*(B*x-offset)
     * to the original primal energy 0.5*x^T*K*x-f^T*x, with optional A*x <= b.
     * The combined dual is solved for lambda >= 0 and lower <= q <= upper;
     * K must be positive definite. B uses original global primal columns, and
     * offset/lower/upper use its distributed row layout. Bounds may be asymmetric
     * or fixed; use PETSC_NINFINITY/PETSC_INFINITY for unbounded endpoints.
     * Retains PETSc references; callers may update values between solves.
     */
    void set_dual_box(Mat B, Vec offset, Vec lower, Vec upper);

    /** Borrow lambda in A row numbering, valid until the next solve or constraint setter.
     * Null before a combined solve or when A has no rows. */
    Vec dual_ineq_multipliers() const { return ineq_multipliers_; }

    /** Borrow q in B row numbering, valid until the next solve or set_dual_box().
     * K*x + A^T*lambda + B^T*q = f after a combined solve. */
    Vec dual_box_multipliers() const { return dual_box_multipliers_; }

    LinSys_PETSC::SolveInfo solve() override;

    /**
     * Returns information on absolute solver accuracy
     */
    inline double get_absolute_accuracy() override {
       return a_tol_;
    };

    void view(string text="") override;

    double get_solution_precision() override;

    /**
     * Compute and return residual:
     * - For the combined dual solve return the original primal problem's algebraic
     *   KKT norm: equilibrium, inequality complementarity, and box projection
     *   residual, with unit projection steps.
     * - The inequality-only path retains its K*x-f residual.
     */
    double compute_residual() override;

    ~LinSys_PERMON( );

private:
    /**
     * Solve dual quadratic problem of the form
         min 0.5*z^T*D*K^{-1}*D^T*z - (D*K^{-1}*f_shifted-c)^T*z,
     * with D=[A; B_active], c=[b; offset_active] in MatNest distributed row order.
     * Fixed box multipliers are eliminated in f_shifted=f-B_fixed^T*q_fixed.
     * Recover x in original global primal numbering and both full multiplier blocks.
     */
    LinSys::SolveInfo solve_dual_qp();

    /// Registrar of class to factory
    static const int registrar;

    // make a pointer to the data array out of a std::vector
    template<typename T>
    T *  makePetscPointer_( std::vector<T> & array )
    {
        if ( array.size() ) return &(array[0]);
        return PETSC_NULLPTR;
    }

    // PetscScalar to double casting functor
    struct PetscScalar2Double_ : public std::unary_function< PetscScalar, double >
    {
        double operator()( PetscScalar arg ) 
        {
            return static_cast<double>( arg );
        }
    };

protected:

    Mat matrix_ineq_;                     ///< Borrowed A in original global primal column indexing.
    Vec ineq_;                            ///< Borrowed b in the distributed A row layout.
    Mat matrix_dual_box_ = nullptr;       ///< Retained B in original global primal column indexing.
    Vec dual_box_offset_ = nullptr;       ///< Retained offset in the distributed B row layout.
    Vec dual_box_lower_ = nullptr;        ///< Retained lower multiplier bounds in the B row layout.
    Vec dual_box_upper_ = nullptr;        ///< Retained upper multiplier bounds in the B row layout.
    Vec ineq_multipliers_ = nullptr;      ///< Owned lambda in the A row layout, replaced on solve/set.
    Vec dual_box_multipliers_ = nullptr;  ///< Owned q in the B row layout, replaced on solve/set.

    QP                 system;
    QPS                solver;
};

#endif /* LA_LINSYS_PERMON_HH_ */
