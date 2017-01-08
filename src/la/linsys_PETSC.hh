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
 * @file    linsys_PETSC.hh
 * @brief   Solver based on the original PETSc solver using MPIAIJ matrix and succesive Schur complement construction
 * @author  Jakub Sistek
 */

#ifndef LA_LINSYS_PETSC_HH_
#define LA_LINSYS_PETSC_HH_

// derived from base linsys
#include "la/linsys.hh"

#include "la/distribution.hh"
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"

#include "petscksp.h"

class LinSys_PETSC : public LinSys
{

public:
	typedef LinSys FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    LinSys_PETSC(const  Distribution * rows_ds);

    /**
     * Copy constructor.
     */
    LinSys_PETSC( LinSys_PETSC &other );


    void set_tolerances(double  r_tol, double a_tol, unsigned int max_it) override;

    /**
     * Returns whole Distribution class for distribution of the solution.
     */
    inline const Distribution* get_ds( )
    { 
        return rows_ds_; 
    }

    const Mat *get_matrix() override
    { 
        return &matrix_;
    }

    const Vec *get_rhs() override
    { 
        return &rhs_;
    }

    PetscErrorCode set_matrix(Mat &matrix, MatStructure str) override
    {
        matrix_changed_ = true;
    	return MatCopy(matrix, matrix_, str);
    }

    PetscErrorCode set_rhs(Vec &rhs) override
    {
        rhs_changed_ = true;
    	return VecCopy(rhs, rhs_);
    }

    PetscErrorCode mat_zero_entries() override
    {
        matrix_changed_ = true;
    	return MatZeroEntries(matrix_);
    }

    PetscErrorCode rhs_zero_entries() override
    {
        rhs_changed_ = true;
    	return VecSet(rhs_, 0);
    }

    void start_allocation() override;

    void start_add_assembly() override;

    void start_insert_assembly() override;

    void mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals ) override;

    void rhs_set_values( int nrow, int *rows, double *vals ) override;

    void preallocate_values(int nrow,int *rows,int ncol,int *cols);

    void preallocate_matrix();

    void finish_assembly() override;

    void finish_assembly( MatAssemblyType assembly_type );

    void apply_constrains( double scalar = 1. ) override;

    void set_initial_guess_nonzero(bool set_nonzero = true);

    int solve() override;

    /**
     * Returns information on absolute solver accuracy
     */
    inline double get_absolute_accuracy() override {
       return a_tol_;
    };

    void view( ) override;

    /**
     * Sets specific parameters of LinSys_PETSC defined by user in input file and used to calculate
     */
    void set_from_input(const Input::Record in_rec) override;

    double get_solution_precision() override;

    double compute_residual() override;

    ~LinSys_PETSC( );

private:
    /// Registrar of class to factory
    static const int registrar;

    // make a pointer to the data array out of a std::vector
    template<typename T> 
    T *  makePetscPointer_( std::vector<T> & array )
    {
        if ( array.size() ) return &(array[0]);
        return PETSC_NULL;
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

    std::string params_;		 //!< command-line-like options for the PETSc solver

    bool    init_guess_nonzero;  //!< flag for starting from nonzero guess

    Mat     matrix_;             //!< Petsc matrix of the problem.
    Vec     rhs_;                //!< PETSc vector constructed with vx array.
    Vec     residual_;

    double  *v_rhs_;             //!< local RHS array pointing to Vec rhs_

    Vec     on_vec_;             //!< Vectors for counting non-zero entries in diagonal block.
    Vec     off_vec_;            //!< Vectors for counting non-zero entries in off-diagonal block.


    double  solution_precision_; // precision of KSP system solver

    KSP                system;
    KSPConvergedReason reason;


};

#endif /* LA_LINSYS_PETSC_HH_ */
