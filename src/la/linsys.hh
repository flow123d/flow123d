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
 * @file
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
#include "input/input_type.hh"
#include "input/accessors.hh"


#include <mpi.h>

#include <vector>

// PETSc includes
#include "petscmat.h"
//#include "petscvec.h"
//#include "petscksp.h"


class LinSys
{
friend class SchurComplement;
public:
    // Abstract Input Record for LinSys initialization
    static Input::Type::AbstractRecord input_type;

    typedef enum {
        INSERT=INSERT_VALUES,
        ADD=ADD_VALUES,
        ALLOCATE,
        DONE,
        NONE
    } SetValuesMode;

    /*typedef enum {
        PETSC,
        BDDC
        //PETSC_schur_complement   // possibly we can implement Schur as another kind of lin solver
        //PETSC_MPIAIJ_preallocate_by_assembly,
        //PETSC_MPIAIJ_assembly_by_triples,
    } LinSysType;*/

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
    LinSys( Distribution * rows_ds,
            MPI_Comm comm = MPI_COMM_WORLD )
      : lsize_( rows_ds->lsize() ), rows_ds_(rows_ds), comm_( comm ), solution_(NULL), v_solution_(NULL),
        positive_definite_( false ), negative_definite_( false ), symmetric_( false ),
        spd_via_symmetric_general_( false ), status_( NONE )
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
      spd_via_symmetric_general_(other.spd_via_symmetric_general_), globalSolution_(other.globalSolution_),
      constraints_(other.constraints_), residual_norm_(other.residual_norm_), in_rec_(other.in_rec_)
    {
    	ASSERT( false, "Using copy constructor of LinSys is not allowed!");
    	set_solution(other.v_solution_);
    };

    // Particular type of the linear system.
    //LinSysType type;  //!< anyone can inquire my type

    virtual void load_mesh( const int nDim, const int numNodes, const int numDofs,
                            const std::vector<int> & inet, 
                            const std::vector<int> & nnet, 
                            const std::vector<int> & nndf, 
                            const std::vector<int> & isegn, 
                            const std::vector<int> & isngn, 
                            const std::vector<int> & isvgvn,
                            const std::vector<double> & xyz,
                            const std::vector<double> & element_permeability,
                            const int meshDim )
    {
        ASSERT( false, "Function load_mesh is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    virtual void load_diagonal( std::map<int,double> & diag )
    {
        ASSERT( false, "Function load_diagonal is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

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
     */
    virtual const Mat &get_matrix()
    {
        ASSERT( false, "Function get_matrix is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Returns RHS vector  (only for PETSC solvers)
     */
    virtual const Vec &get_rhs()
    {
        ASSERT( false, "Function get_rhs is not implemented for linsys type %s \n.", typeid(*this).name() );
    }
    
    /**
     * Sets matrix changed flag  (only for PETSC solvers)
     */
    virtual void set_matrix_changed()
    {
        ASSERT( false, "Function set_matrix_changed is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Sets rhs changed flag  (only for PETSC solvers)
     */
    virtual void set_rhs_changed()
    {
        ASSERT( false, "Function set_rhs_changed is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Sets PETSC matrix (only for PETSC solvers)
     */
    virtual PetscErrorCode set_matrix(Mat &matrix, MatStructure str)
    {
        ASSERT( false, "Function set_matrix is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Sets RHS vector  (only for PETSC solvers)
     */
    virtual PetscErrorCode set_rhs(Vec &rhs)
    {
        ASSERT( false, "Function set_rhs is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    virtual PetscErrorCode mat_zero_entries()
    {
    	ASSERT( false, "Function mat_zero_entries is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    virtual PetscErrorCode rhs_zero_entries()
    {
    	ASSERT( false, "Function vec_zero_entries is not implemented for linsys type %s \n.", typeid(*this).name() );
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
     * Returns whole solution vector.
     */
    virtual void get_whole_solution( std::vector<double> & globalSolution )
    {
        ASSERT( false, "Function get_whole_solution is not implemented for linsys type %s \n.", typeid(*this).name() );
    }

    /**
     * Inserts solution vector.
     */
    virtual void set_whole_solution( std::vector<double> & globalSolution )
    {
        ASSERT( false, "Function set_whole_solution is not implemented for linsys type %s \n.", typeid(*this).name() );
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

    /**
     * Shortcut to assembly into matrix and RHS in one call.
     * This can also apply constrains at assembly time (only in add assembly regime).
     *
     * Constrains can either be set before through add_constraint. Or by additional parameters if we
     * have only per element knowledge about boundary conditions.
     *
     */
    void set_values( int nrow,int *rows,int ncol,int *cols,double *mat_vals, double *rhs_vals )
//                            std::vector<bool> &constrains_row_mask=std::vector<bool>(), double * constrain_values=NULL )
    {
        mat_set_values(nrow, rows, ncol, cols, mat_vals);
        rhs_set_values(nrow, rows, rhs_vals);
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
     * Returns information on relative solver accuracy
     */
    double get_relative_accuracy(){
       return r_tol_;
    };

    /**
     * Returns information on absolute solver accuracy
     */
    virtual double get_absolute_accuracy(){
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
    //virtual void view(std::ostream output_stream, int * output_mapping = NULL)
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

    ~LinSys()
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

    Vec      solution_;          //!< PETSc vector constructed with vb array.
    double  *v_solution_;        //!< local solution array pointing into Vec solution_
    bool     own_solution_;      //!< Indicates if the solution array has been allocated by this class

    double  residual_norm_;      //!< local solution array pointing into Vec solution_

    ConstraintVec_   constraints_;

    std::vector<double>  globalSolution_; //!< global solution in numbering for linear system

    Input::Record in_rec_;        // structure contained parameters of LinSys defined in input file

};

#endif /* LA_LINSYS_HH_ */
