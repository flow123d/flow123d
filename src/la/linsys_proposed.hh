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

#include "petscmat.h"
#include "private/matimpl.h"

#include "la/distribution.hh"


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
 */

class LinSys
{
public:
    typedef enum {
        INSERT=INSERT_VALUES,
        ADD=ADD_VALUES,
        ALLOCATE,
        DONE,
        NONE
    } SetValuesMode;

    typedef enum {
        PETSC_MPIAIJ_preallocate_by_assembly,
        PETSC_MPIAIJ_assembly_by_triples,
        BDDC,
        PETSC_schur_complement   // possibly we can implement Schur as another kind of lin solver
    } LinSysType;

    /**
     * Constructor.
     * Constructor of abstract class should not be called directly, but is used for initialization of member common
     * to all particular solvers.
     *
     * @param lsize - local size of the solution vector
     * @param sol_array - optionally one can provide pointer to array allocated to size lsize, where
     *  the solution should be stored,
     */
    LinSys(unsigned int lsize, double *sol_array = NULL);

    /// Particular type of the linear system.
    LinSysType  type;   ///< MAT_IS or MAT_MPIAIJ anyone can inquire my type


    /**
     *  Returns global system size.
     */
    inline unsigned int size()
    { return vec_ds.size(); }

    /**
     * Returns local system size. (for partitioning of solution vectors)
     * for PETSC_MPIAIJ it is also partitioning of the matrix
     */
    inline unsigned int vec_lsize()
    { return vec_ds.lsize(); }

    /**
     * Returns whole Distribution class for distribution of the solution.
     */
    inline const Distribution &ds()
    { return vec_ds; }

    /**
     * Returns PETSC matrix (only for PETSC solvers)
     */
    inline const Mat &get_matrix()
    { return matrix; }

    /**
     * Returns RHS vector  (only for PETSC solvers)
     */
    inline const Vec &get_rhs()
    { return rhs; }

    /**
     *  Returns PETSC vector with solution. Underlying array can be provided on construction.
     *  Can this be implemented for BDDC?
     */
    inline const Vec &get_solution()
    { return solution; }

    /**
     * Returns local part of solution vector.
     * Can this be implemented for BDDC?
     */
    inline double *get_solution_array()
    { return v_solution; }

    /**
     * Switch linear system into allocating assembly. (only for PETSC_MPIAIJ_preallocate_by_assembly)
     */
    virtual void start_allocation()=0;

    /**
     * Switch linear system into adding assembly. (the only one supported by triplets ??)
     */
    void start_add_assembly();
    /**
     * Switch linear system into insert assembly. (not currently used)
     */
    void start_insert_assembly();


    /**
     * May not be necessary. For PETSC this should call MatEndAssembly with MAT_PARTIAL_ASSEMBLY
     */
    void partial_assembly();
    /**
     * Finish assembly of the whole system. For PETSC this should call MatEndAssembly with MAT_FINAL_ASSEMBLY
     */
    void finish_assembly();

    /**
     *  Assembly full rectangular submatrix into the system matrix.
     *  Should be virtual, implemented differently in  particular solvers.
     */
    virtual void mat_set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *vals) = 0;


    /**
     * Shortcut for assembling just one element into the matrix.
     * Similarly we can provide method accepting armadillo matrices.
     */
    inline void mat_set_value(int row,int col,PetscScalar val)
    { mat_set_values(1,&row,1,&col,&val); }

    /**
     *  Set values of the system right-hand side.
     *  Should be virtual, implemented differently in  particular solvers.
     */
    inline void rhs_set_values(int nrow,int *rows,PetscScalar *vals)
    {
        switch (status) {
            case INSERT:
            case ADD:
                VecSetValues(rhs,nrow,rows,vals,(InsertMode)status); break;
            case ALLOCATE: break;
            default: ASSERT(0, "LinSys's status disallow set values.\n");
        }
    }

    /**
     * Shorcut for assembling just one element into RHS vector.
     */
    inline void rhs_set_value(int row,PetscScalar val)
    { rhs_set_values(1,&row,&val); }

    /**
     * Shortcut to assembly into matrix and RHS in one call.
     * This can also apply constrains at assembly time (only in add assembly regime).
     *
     * Constrains can either be set before through add_constrain. Or by additional parameters if we
     * have only per element knowledge about boundary conditions.
     *
     */
    inline void set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *mat_vals, PetscScalar *rhs_vals,
            std::vector<bool> &constrains_row_mask=std::vector(0), double * constrain_values=NULL)
    {
        mat_set_values(nrow, rows, ncol, cols, mat_vals);
        rhs_set_values(nrow, rows, rhs_vals);
    }

    /**
     * Adds Dirichlet constrain.
     * @param row - global numeb of row that should be eliminated.
     * @param value - solution value at the given row
     */
    void add_constrain(int row, double value);

    /**
     * Apply constrains to assembled matrix. Constrains are given by pairs: global row index, value.
     * i.e. typedef pair<unsigned int, double> Constrain;
     *
     * What is th meaning of ( const double factor ) form Cambridge code?
     */

    void apply_constrains(std::vector<Constrain> & constraints);


    /**
     * Solve the system.
     *
     * parameters should by provided in input file (currently INI file, but will be changed to JSON)
     *
     * If we realize that we need to set something, rather add some set_* function.
     *
     * double tol = 1.e-7,                        //!< tolerance on relative residual ||res||/||rhs||
       int  numLevels = 2,                        //!< number of levels
       std::vector<int> *  numSubAtLevels = NULL, //!< number of subdomains at levels
       int verboseLevel = 0,                      //!< level of verbosity of BDDCML library
                                                  //!< ( 0 - only fatal errors reported,
                                                  //!<   1 - mild output,
                                                  //!<   2 - detailed output )
       int  maxIt = 1000,                         //!< maximum number of iterations
       int  ndecrMax = 30 );                      //!< maximum number of iterations with non-decreasing residual
                                                  //!< ( used to stop diverging process )
     *
     *
     * Returns convergence reason (form PETSC)
     */
    virtual int solve();


    /**
     * Provides user knowledge about symmetry.
     */
    inline void set_symmetric(bool flag = true)
    {
        symmetric = flag;
        if (!flag) set_positive_definite(false);
    }

    inline bool is_symmetric()
    { return symmetric; }

    /**
     * Provides user knowledge about positive definiteness.
     */
    inline void set_positive_definite(bool flag = true)
    {
        positive_definite = flag;
        if (flag) set_symmetric();
    }

    inline bool is_positive_definite()
    { return positive_definite; }


    /**
     *  Output the system in the Matlab format possibly with given ordering.
     *  Rather we shoud provide output operator <<, since it is more flexible.
     */
    void view(std::ostream output_stream, int * output_mapping = NULL);

    virtual ~LinSys();

protected:
    /**
     * Protected methods used in preallocate_by_assembly solvers.
     */
    virtual void preallocate_matrix()=0;
    /**
     * Protected methods used in preallocate_by_assembly solvers.
     */
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols)=0;


    Distribution vec_ds;        ///< Distribution of continuous blocks of system rows among the processors.
    bool     symmetric;         ///< Flag for the symmetric system.
    bool     positive_definite; ///< Flag for positive definite system.
    bool     own_solution;      ///< Indicates if the solution array has been allocated by this class.
    SetValuesMode status;       ///< Set value status of the linear system.

    Mat     matrix;             ///< Petsc matrix of the problem.
    Vec     rhs;                ///< PETSc vector constructed with vx array.
    Vec     solution;           ///< PETSc vector constructed with vb array.
    double  *v_rhs;             ///< RHS vector.
    double  *v_solution;        ///< Vector of solution.

    // for MATIS
    int *subdomain_indices;     ///< Remember indices which created mapping
    Mat local_matrix;           ///< local matrix of subdomain (used in LinSys_MATIS)

    friend void SchurComplement::form_schur();
    friend class SchurComplement;
};



class LinSys_MPIAIJ : public LinSys
{
public:
    LinSys_MPIAIJ(unsigned int lsize, double *sol_array=NULL) : LinSys(lsize, sol_array) {};
    virtual void start_allocation();
    virtual void preallocate_matrix();
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols);
    virtual void view_local_matrix();
    virtual ~LinSys_MPIAIJ();

private:

    Vec     on_vec,off_vec; ///< Vectors for counting non-zero entries.

};


class LinSys_MATIS : public LinSys
{
public:
    LinSys_MATIS(unsigned int lsize, int sz, int *global_row_4_sub_row, double *sol_array=NULL);
    virtual void start_allocation();
    virtual void preallocate_matrix();
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols);
    virtual void view_local_matrix();

    inline VecScatter get_scatter()
    { return sub_scatter; }
    /// Get local subdomain size.
    inline int get_subdomain_size()
    { return subdomain_size; }
  
    virtual ~LinSys_MATIS();

private:

    ISLocalToGlobalMapping map_local_to_global; ///< PETSC mapping form local indexes of subdomain to global indexes
    int loc_rows_size;                          ///<
    int *loc_rows;                              ///< Small auxiliary array for translation of global indexes to local
                                                ///< during preallocate_set_values. However for MatSetValues
    VecScatter sub_scatter;                     ///< from global vector with no overlaps constructs local (subdomain)
                                                ///< vectors with overlaps
                                                ///< copy of scatter created and used in PETSc in Mat_IS

    int subdomain_size;                         ///< size of subdomain in MATIS matrix
    int *subdomain_nz;                          ///< For counting non-zero enteries of local subdomain.

    // mimic PETSc struct for IS matrices - included from matis.h
    // used to access private PETSc data
    typedef struct {
       Mat                    A;             /* the local Neumann matrix */
       VecScatter             ctx;           /* update ghost points for matrix vector product */
       Vec                    x,y;           /* work space for ghost values for matrix vector product */
       ISLocalToGlobalMapping mapping;
       int                    rstart,rend;   /* local row ownership */
       PetscTruth             pure_neumann;
    } MatMyIS ;
};



#endif /* LA_LINSYS_HH_ */
