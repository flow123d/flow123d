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
 * $Id: la_linsys.hh 1299 2011-08-23 21:42:50Z jakub.sistek $
 * $Revision: 1299 $
 * $LastChangedBy: jakub.sistek $
 * $LastChangedDate: 2011-08-23 23:42:50 +0200 (Tue, 23 Aug 2011) $
 *
 * @file
 * @brief   Solver based on the original PETSc solver using MPIAIJ matrix and succesive Schur complement construction
 * @author  Jakub Sistek
 *
 *
 */

#ifndef LA_LINSYS_PETSC_HH_
#define LA_LINSYS_PETSC_HH_

// derived from base linsys
#include "la/linsys.hh"

#include "la/distribution.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"

class LinSys_PETSC : public LinSys
{

public:
    static Input::Type::Record input_type;

    LinSys_PETSC( Distribution * rows_ds,
                  const MPI_Comm comm = PETSC_COMM_WORLD ); 

    /**
     * Returns whole Distribution class for distribution of the solution.
     */
    inline const Distribution* get_ds( )
    { 
        return rows_ds_; 
    }

    const Mat &get_matrix()
    { 
        return matrix_; 
    }

    const Vec &get_rhs()
    { 
        return rhs_; 
    }

    PetscErrorCode set_matrix(Mat &matrix, MatStructure str)
    {
    	return MatCopy(matrix, matrix_, str);
    }

    PetscErrorCode set_rhs(Vec &rhs)
    {
    	return VecCopy(rhs, rhs_);
    }

    PetscErrorCode mat_zero_entries()
    {
    	return MatZeroEntries(matrix_);
    }

    PetscErrorCode rhs_zero_entries()
    {
    	return VecSet(rhs_, 0);
    }

    void start_allocation();

    void start_add_assembly();

    void start_insert_assembly();

    void mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals );

    void rhs_set_values( int nrow, int *rows, double *vals );

    void preallocate_values(int nrow,int *rows,int ncol,int *cols);

    void preallocate_matrix();

    void finish_assembly();

    void finish_assembly( MatAssemblyType assembly_type );

    void apply_constrains( double scalar = 1. );

    void set_initial_guess_nonzero(bool set_nonzero = true);

    int solve();

    /**
     * Returns information on absolute solver accuracy
     */
    inline double get_absolute_accuracy(){
       return a_tol_;
    };

    void get_whole_solution( std::vector<double> & globalSolution );

    void view( );

    /**
     * Sets specific parameters of LinSys_PETSC defined by user in input file and used to calculate
     */
    void set_from_input(const Input::Record in_rec);


    ~LinSys_PETSC( );

private:
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

    void gatherSolution_( );

private:

    std::string params_;		 //!< command-line-like options for the PETSc solver

    bool    init_guess_nonzero;  //!< flag for starting from nonzero guess

    Mat     matrix_;             //!< Petsc matrix of the problem.
    Vec     rhs_;                //!< PETSc vector constructed with vx array.

    double  *v_rhs_;             //!< local RHS array pointing to Vec rhs_

    Vec     on_vec_;             //!< Vectors for counting non-zero entries in diagonal block.
    Vec     off_vec_;            //!< Vectors for counting non-zero entries in off-diagonal block.

};

#endif /* LA_LINSYS_PETSC_HH_ */
