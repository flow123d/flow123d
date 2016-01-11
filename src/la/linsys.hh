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
 * @file    linsys.hh
 * @brief   Wrappers for linear systems based on MPIAIJ and MATIS format.
 * @author  Jan Brezina
 *
 * Linear system only keep together matrix, RHS and the solution vector.
 *
 */

//=============================================================================
//
// LINER SYSTEM ROUTINES - linear system use wrapers of PETSc assembling routins
// in order to allow counting of allocation and filling of matrix to be done by the same code
//
// we want to allow allocations of nonlocal rows, to this end we count on- and off-processor
// members in an parallel Vector
//
//=============================================================================

#ifndef LA_LINSYS_HH_
#define LA_LINSYS_HH_

/**
 *  @brief  Abstract linear system class.
 *
 *  Linear system consists of Matrix, RHS and solution.
 *  It provides general methods for:
 *  - matrix preallocation
 *  - assembling matrix and RHS
 *  - application of linear constrains (Dirichlet conditions) either during assembly or
 *    on assembled system
 *  - solution of the system
 *  - output in Matlab format
 *
 *  Method operates on the system as single object. But some methods for matrix only manipulation
 *  can be provided until we have matrix as separate class.
 *
 *  TODO:
 *  - simplify constructors
 *  - introduce set_solution_vector() and set_rhs() in order to solve multiple systems with same matrix
 *  - simplify constructor, make common interface, rather introduce particular seters for particular solvers
 *  - which parameters expose through Input::Record
 *  - why we need lsize in constructors if we have Distribution
 */

#include "system/global_defs.h"
#include "la/distribution.hh"
#include "input/input_type_forward.hh"
#include "input/accessors.hh"


#include <mpi.h>

#include <vector>
#include <armadillo>

// PETSc includes
#include "petscmat.h"


class LinSys
{
friend class SchurComplement;
public:
    // Abstract Input Record for LinSys initialization
    static Input::Type::Abstract & get_input_type();

    typedef enum {
        INSERT=INSERT_VALUES,
        ADD=ADD_VALUES,
        ALLOCATE,
        DONE,
        NONE
    } SetValuesMode;

protected:
    typedef std::pair<unsigned,double>       Constraint_;
    typedef std::vector< Constraint_ >       ConstraintVec_;

public:
    /**
     * Constructor.
     * Constructor of abstract class should not be called directly, but is used for initialization of member common
     * to all particular solvers.
     *
     * @param comm - MPI communicator
     *
     * TODO: Vector solution_ is now initialized to NULL, but it should be rather allocated
     * in the constructor instead of the method set_solution().
     */
    LinSys(const  Distribution *rows_ds)
      : comm_( rows_ds->get_comm() ), status_( NONE ), lsize_( rows_ds->lsize() ), rows_ds_(rows_ds),
        symmetric_( false ), positive_definite_( false ), negative_definite_( false ),
        spd_via_symmetric_general_( false ), solution_(NULL), v_solution_(NULL)
    { 
        int lsizeInt = static_cast<int>( rows_ds->lsize() );
        int sizeInt;
        MPI_Allreduce ( &lsizeInt, &sizeInt, 1, MPI_INT, MPI_SUM, comm_ );
        size_ = static_cast<unsigned>( sizeInt );

    };

    /**
     * Copy constructor.
     */
    LinSys(LinSys &other)
    : r_tol_(other.r_tol_), a_tol_(other.a_tol_), max_it_(other.max_it_), comm_(other.comm_), status_(other.status_),
      lsize_( other.rows_ds_->lsize() ), size_(other.size_), rows_ds_(other.rows_ds_), symmetric_(other.symmetric_),
      positive_definite_(other.positive_definite_), negative_definite_( other.negative_definite_ ),
      spd_via_symmetric_general_(other.spd_via_symmetric_general_), matrix_changed_(other.matrix_changed_),
	  rhs_changed_(other.rhs_changed_), residual_norm_(other.residual_norm_), constraints_(other.constraints_),
      globalSolution_(other.globalSolution_), in_rec_(other.in_rec_)

    {
    	ASSERT( false, "Using copy constructor of LinSys is not allowed!");
    	set_solution(other.v_solution_);
    };

    /**
     *  Returns global system size.
     */
    inline unsigned int size()
    { 
        return size_; 
    }

    /**
     * Returns local system size. (for partitioning of solution vectors)
     * for PETSC_MPIAIJ it is also partitioning of the matrix
     */
    inline unsigned int vec_lsize()
    { 
        return lsize_; 
    }

    /**
     * Returns PETSC matrix (only for PETSC solvers)
     *
     * If matrix is changed, method set_matrix_changed() must be called.
     * Example:
	 * @CODE
	 *   MatDiagonalSet(schur->get_matrix(), new_diagonal, ADD_VALUES);
	 *   schur->set_matrix_changed();
	 * @ENDCODE
     */
    virtual const Mat *get_matrix()
    {
        ASSERT( false, "Function get_matrix is not implemented for linsys type %s \n.", typeid(*this).name() );
        return NULL;
    }

    /**
     * Returns RHS vector  (only for PETSC solvers)
     *
     * If vector is changed, method set_rhs_changed() must be called.
     * Example:
	 * @CODE
	 *   VecScale(schur->get_rhs(), -1.0);
	 *   schur->set_rhs_changed();
	 * @ENDCODE
     */
    virtual const Vec *get_rhs()
    {
        ASSERT( false, "Function get_rhs is not implemented for linsys type %s \n.", typeid(*this).name() );
        return NULL;
    }
    
    /**
     * Sets matrix changed flag.
     */
    void set_matrix_changed()
    { matrix_changed_ = true;}

    /**
     * Sets rhs changed flag  (only for PETSC solvers)
     */
    void set_rhs_changed()
    { rhs_changed_ = true; }

    /**
     * Returns true if the system matrix has changed since the last solve.
     */
    bool is_matrix_changed()
    { return matrix_changed_;}

    /**
     * Returns true if the system RHS has changed since the last solve.
     */
    bool is_rhs_changed()
    { return rhs_changed_;}


    /**
     * Sets PETSC matrix (only for PETSC solvers)
     */
    virtual PetscErrorCode set_matrix(Mat &matrix, MatStructure str)
    {
        ASSERT( false, "Function set_matrix is not implemented for linsys type %s \n.", typeid(*this).name() );
        return 0;
    }

    /**
     * Sets RHS vector  (only for PETSC solvers)
     */
    virtual PetscErrorCode set_rhs(Vec &rhs)
    {
        ASSERT( false, "Function set_rhs is not implemented for linsys type %s \n.", typeid(*this).name() );
        return 0;
    }

    /**
     * Clears entries of the matrix 
     */
    virtual PetscErrorCode mat_zero_entries()
    {
    	ASSERT( false, "Function mat_zero_entries is not implemented for linsys type %s \n.", typeid(*this).name() );
    	return 0;
    }

    /**
     * Clears entries of the right-hand side 
     */
    virtual PetscErrorCode rhs_zero_entries()
    {
    	ASSERT( false, "Function vec_zero_entries is not implemented for linsys type %s \n.", typeid(*this).name() );
    	return 0;
    }

    /**
     *  Returns PETSC vector with solution. Underlying array can be provided on construction.
     */
    const Vec &get_solution()
    { 
        return solution_; 
    }

    /**
     * Create PETSc solution
     */
    void set_solution(double *sol_array) {
        if (sol_array == NULL) {
            v_solution_   = new double[ rows_ds_->lsize() + 1 ];
            own_solution_ = true;
        }
        else {
            v_solution_ = sol_array;
            own_solution_ = false;
        }
        PetscErrorCode ierr;
        ierr = VecCreateMPIWithArray( comm_,1, rows_ds_->lsize(), PETSC_DECIDE, v_solution_, &solution_ ); CHKERRV( ierr );
    }

    /**
     *  Returns PETSC subarray with solution. Underlying array can be provided on construction.
     */
    double *get_solution_array()
    { 
        return v_solution_; 
    }

    /**
     * Switch linear system into allocating assembly. (only for PETSC_MPIAIJ_preallocate_by_assembly)
     */
    virtual void start_allocation()
    {
        ASSERT( false, "Function start_allocation is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Switch linear system into adding assembly. (the only one supported by triplets ??)
     */
    virtual void start_add_assembly()
    {
        ASSERT( false, "Function start_add_assembly is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Switch linear system into insert assembly. (not currently used)
     */
    virtual void start_insert_assembly()
    {
        ASSERT( false, "Function start_insert_assembly is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Finish assembly of the whole system. For PETSC this should call MatEndAssembly with MAT_FINAL_ASSEMBLY
     */
    virtual void finish_assembly( )=0;

    /**
     *  Assembly full rectangular submatrix into the system matrix.
     *  Should be virtual, implemented differently in  particular solvers.
     */
    virtual void mat_set_values(int nrow,int *rows,int ncol,int *cols,double *vals)=0;

    /**
     * Shortcut for assembling just one element into the matrix.
     * Similarly we can provide method accepting armadillo matrices.
     */
    void mat_set_value(int row,int col,double val)
    { mat_set_values(1,&row,1,&col,&val); }

    /**
     *  Set values of the system right-hand side.
     *  Should be virtual, implemented differently in  particular solvers.
     */
    virtual void rhs_set_values(int nrow,int *rows,double *vals)=0;

    /**
     * Shorcut for assembling just one element into RHS vector.
     */
    void rhs_set_value(int row,double val)
    { rhs_set_values(1,&row,&val); }


    /// Set values in the system matrix and values in the right-hand side vector on corresponding rows.
    inline void set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *mat_vals, PetscScalar *rhs_vals)
    {
        mat_set_values(nrow, rows, ncol, cols, mat_vals);
        rhs_set_values(nrow, rows, rhs_vals);
    }

    /**
     * Shortcut to assembly into matrix and RHS in one call, possibly apply Dirichlet boundary conditions.
     * @p row_dofs - are global indices of rows of dense @p matrix and rows of dense vector @rhs in global system
     * @p col_dofs - are global indices of columns of the matrix, and possibly
     *
     * Application of Dirichlet conditions:
     * 1) Rows with negative dofs are set to zero.
     * 2) Cols with negative dofs are eliminated.
     * 3) If there are entries on global diagonal. We determine value K either from diagonal of local matrix, or (if it is zero) from
     *    diagonal average.
     *
     * Caveats:
     * - can not set dirichlet condition on zero dof 
     * - Armadillo stores matrix in column first form (Fortran like) which makes it not well suited 
     *   for passing local matrices.
     *
     */
    void set_values(std::vector<int> &row_dofs, std::vector<int> &col_dofs,
    		        const arma::mat &matrix, const arma::vec &rhs,
    		        const arma::vec &row_solution, const arma::vec &col_solution)

    {
    	arma::mat tmp = matrix.t();
    	arma::vec tmp_rhs = rhs;
    	bool negative_row = false;
    	bool negative_col = false;

    	for(unsigned int l_row = 0; l_row < row_dofs.size(); l_row++)
    		if (row_dofs[l_row] < 0) {
                        tmp_rhs(l_row)=0.0;
    			tmp.col(l_row).zeros();
    			negative_row=true;
    		}

    	for(unsigned int l_col = 0; l_col < col_dofs.size(); l_col++)
    		if (col_dofs[l_col] < 0) {
    			tmp_rhs -= matrix.col(l_col) * col_solution[l_col];
    			tmp.row(l_col).zeros();
                        negative_col=true;
    		}
    		

    	if (negative_row && negative_col) {
    		// look for diagonal entry
        	for(unsigned int l_row = 0; l_row < row_dofs.size(); l_row++)
        		if (row_dofs[l_row] < 0)
        	    	for(unsigned int l_col = 0; l_col < col_dofs.size(); l_col++)
        	    		if (col_dofs[l_col] < 0 && row_dofs[l_row] == col_dofs[l_col]) {
        	    			double new_diagonal = fabs(matrix.at(l_row,l_col));
        	    			if (new_diagonal == 0.0) {
        	    				if (matrix.is_square()) {
        	    					new_diagonal = arma::sum( abs(matrix.diag())) / matrix.n_rows;
        	    				} else {
        	    					new_diagonal = arma::accu( abs(matrix) ) / matrix.n_elem;
        	    				}
        	    			}
        	    			tmp.at(l_col, l_row) = new_diagonal;
        	    			tmp_rhs(l_row) = new_diagonal * row_solution[l_row];
        	    		}

    	}

    	if (negative_row)
    		for(int &row : row_dofs) row=abs(row);

    	if (negative_col)
    		for(int &col : col_dofs) col=abs(col);


        mat_set_values(row_dofs.size(), const_cast<int *>(&(row_dofs[0])),
        		       col_dofs.size(), const_cast<int *>(&(col_dofs[0])), tmp.memptr() );
        rhs_set_values(row_dofs.size(), const_cast<int *>(&(row_dofs[0])), tmp_rhs.memptr() );
    }

    /**
     * Adds Dirichlet constrain.
     * @param row - global number of row that should be eliminated.
     * @param value - solution value at the given row
     */
    void add_constraint(int row, double value) {

        constraints_.push_back( Constraint_( static_cast<unsigned>( row ), value ) );
    }

    /**
     * Apply constrains to assembled matrix. Constrains are given by pairs: global row index, value.
     * i.e. typedef pair<unsigned int, double> Constrain;
     *
     * What is th meaning of ( const double factor ) form Cambridge code?
     */
    virtual void apply_constrains( double scalar )=0;

    /**
     * Solve the system and return convergence reason.
     */
    virtual int solve()=0;

    /**
     * Returns norm of the initial right hand side
     */
    double get_residual_norm(){
       return residual_norm_;
    };

    /**
     * Explicitly compute residual and its norm for current solution.
     */
    virtual double compute_residual() =0;

    /**
     * Returns information on relative solver accuracy
     */
    double get_relative_accuracy(){
       return r_tol_;
    };

    /**
     * Returns information on absolute solver accuracy
     */
    virtual double get_absolute_accuracy(){
    	return 0.0;
    };

    /**
     * Provides user knowledge about symmetry.
     */
    inline void set_symmetric(bool flag = true)
    {
        symmetric_ = flag;
        if (!flag) {
        	set_positive_definite(false);
        	set_negative_definite(false);
        }
    }

    inline bool is_symmetric()
    { return symmetric_; }

    /**
     * Provides user knowledge about positive definiteness.
     */
    inline void set_positive_definite(bool flag = true)
    {
        positive_definite_ = flag;
        if (flag) {
        	set_symmetric();
        	set_negative_definite(false);
        }
    }

    /**
     * Provides user knowledge about negative definiteness.
     */
    inline void set_negative_definite(bool flag = true)
    {
    	negative_definite_ = flag;
        if (flag) {
        	set_symmetric();
        	set_positive_definite(false);
        }
    }

    inline bool is_positive_definite()
    { return positive_definite_; }

    inline bool is_negative_definite()
    { return negative_definite_; }

    /// TODO: In fact we want to know if the matrix is already preallocated
    /// However to do this we need explicit finalisation of preallocating cycle.
    inline bool is_new() {
        return ( status_ == NONE );
    }

    inline bool is_preallocated() {
        return ( status_ == INSERT || status_ == ADD);
    }

    /**
     * Provides user knowledge about positive definiteness via symmetric general approach.
     * This is useful for solving Darcy flow by mixed hybrid method, where blocks on subdomains are saddle point but 
     * interface among subdomains is only at the block of Lagrange multipliers and is symmetric positive definite.
     * Problem condensed to interface can thus be solved by PCG method, although original problem is saddle point.
     */
    inline void set_spd_via_symmetric_general(bool flag = true)
    {
        spd_via_symmetric_general_ = flag;
        if (flag) set_symmetric();
    }

    inline bool is_spd_via_symmetric_general()
    { return spd_via_symmetric_general_; }


    /**
     *  Output the system in the Matlab format possibly with given ordering.
     *  Rather we shoud provide output operator <<, since it is more flexible.
     */
    virtual void view()
    {
        ASSERT( false, "Function view is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Sets basic parameters of LinSys defined by user in input file and used to calculate
     */
    virtual void set_from_input(const Input::Record in_rec)
    {
    	if (! in_rec.is_empty()) {
    		in_rec_ = Input::Record( in_rec );
    		r_tol_  = in_rec.val<double>("r_tol");
    		max_it_ = in_rec.val<int>("max_it");
    		a_tol_  = 0.01 * r_tol_;
    	}
    }

    /**
     * Get precision of solving
     */
    virtual double get_solution_precision() = 0;

    virtual ~LinSys()
    { 
       PetscErrorCode ierr;
       if ( solution_ ) { ierr = VecDestroy(&solution_); CHKERRV( ierr ); }
       if ( own_solution_ ) delete[] v_solution_;
    }

protected:
    double           r_tol_;  // relative tolerance of linear solver
    double      	 a_tol_;  // absolute tolerance of linear solver
    int              max_it_; // maximum number of iterations of linear solver

    MPI_Comm         comm_;
    SetValuesMode    status_;         //!< Set value status of the linear system.

    const unsigned   lsize_;          //!< local number of matrix rows (non-overlapping division of rows)
    unsigned          size_;          //!< global number of matrix rows, i.e. problem size

    const Distribution * rows_ds_;   //!< final distribution of rows of MH matrix

    bool             symmetric_;
    bool             positive_definite_;
    bool             negative_definite_;
    bool             spd_via_symmetric_general_;

    bool    matrix_changed_;     //!< true if the matrix was changed since the last solve
    bool    rhs_changed_;        //!< true if the right hand side was changed since the last solve

    Vec      solution_;          //!< PETSc vector constructed with vb array.
    double  *v_solution_;        //!< local solution array pointing into Vec solution_
    bool     own_solution_;      //!< Indicates if the solution array has been allocated by this class

    double  residual_norm_;      //!< local solution array pointing into Vec solution_

    ConstraintVec_   constraints_;

    std::vector<double>  globalSolution_; //!< global solution in numbering for linear system

    Input::Record in_rec_;        // structure contained parameters of LinSys defined in input file

};

#endif /* LA_LINSYS_HH_ */
