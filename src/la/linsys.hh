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

#include "petscmat.h"
//#include "private/matimpl.h"

#include "la/schur.hh"
#include "la/distribution.hh"
#include "la/local_to_global_map.hh"


// **************************************************************
/*!  @brief  Linear System structure accepted by Solver module
 *
 *  The system is based on PETSc matrix A and vector of RHS (b) and solution (x),
 *  both vectors are build on the regular arrays vx, vb.
 *  CSR storage are optional and generated on demand by LinSysSetCSR.
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
        MAT_MPIAIJ,
        MAT_IS
    } SetType;

    /// Construct a parallel system with given local size.
    LinSys(unsigned int lsize, double *sol_array = NULL);

    SetType  type;   ///< MAT_IS or MAT_MPIAIJ anyone can inquire my type

    /// @name access members @{
    /// Get global system size.
    inline unsigned int size()
    { return vec_ds.size(); }
    /// Get local system size.
    inline unsigned int vec_lsize()
    { return vec_ds.lsize(); }
    /// Get distribution of rows.
    inline const Distribution &ds()
    { return vec_ds; }
    /// Get matrix.
    inline const Mat &get_matrix()
    { return matrix; }

    /// Get subdomain matrix.
    /*
    inline const Mat &get_matrix_sub()
    { 
       if      (type == MAT_IS)
       {
	  return local_matrix;
       }
    }*/

    /// Get RHS.
    inline const Vec &get_rhs()
    { return rhs; }
    /// Get solution.
    inline const Vec &get_solution()
    { return solution; }
    inline double *get_solution_array()
    { return v_solution; }
    /// @}

    virtual void start_allocation()=0;
    void start_add_assembly();
    void start_insert_assembly();
    void finalize(MatAssemblyType assembly_type=MAT_FINAL_ASSEMBLY);

    virtual void preallocate_matrix()=0;
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols)=0;

    virtual void view_local_matrix()=0;

    /// Set full rectangular submatrix of the system matrix.
    inline void mat_set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *vals)
    {
        switch (status) {
            case INSERT:
            case ADD:
                MatSetValues(matrix,nrow,rows,ncol,cols,vals,(InsertMode)status); break;
            case ALLOCATE:
                preallocate_values(nrow,rows,ncol,cols); break;
            default: DBGMSG("LS SetValues with non allowed insert mode.\n");
        }
    }

    /// Set one element of the system matrix.

    inline void mat_set_value(int row,int col,PetscScalar val)
    { mat_set_values(1,&row,1,&col,&val); }

    /// Set values of the system right-hand side.
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

    /// Set one value in the right-hand side.
    inline void rhs_set_value(int row,PetscScalar val)
    { rhs_set_values(1,&row,&val); }

    /// Set values in the system matrix and values in the right-hand side vector on corresponding rows.
    inline void set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *mat_vals, PetscScalar *rhs_vals)
    {
        mat_set_values(nrow, rows, ncol, cols, mat_vals);
        rhs_set_values(nrow, rows, rhs_vals);
    }

    inline void set_symmetric(bool flag = true)
    {
        symmetric = flag;
        if (!flag) set_positive_definite(false);
    }


    inline bool is_symmetric()
    { return symmetric; }

    inline void set_positive_definite(bool flag = true)
    {
        positive_definite = flag;
        if (flag) set_symmetric();
    }

    inline bool is_positive_definite()
    { return positive_definite; }

    /// TODO: In fact we want to know if the matrix is already preallocated
    /// However to do this we need explicit finalisation of preallocating cycle.
    inline bool is_new() {
        return ( status == NONE );
    }

    inline bool is_preallocated() {
        return (status == INSERT || status == ADD);
    }

    /// Output the system in the Matlab format possibly with given ordering.
    void view(std::ostream output_stream, int * output_mapping = NULL);

    virtual ~LinSys();

protected:
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
    LinSys_MATIS(boost::shared_ptr<LocalToGlobalMap> global_row_4_sub_row, double *sol_array=NULL);
    virtual void start_allocation();
    virtual void preallocate_matrix();
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols);
    virtual void view_local_matrix();

    inline VecScatter get_scatter()
    { return sub_scatter; }
    /// Get local subdomain size.
    inline int get_subdomain_size()
    { return lg_map->size(); }
  
    virtual ~LinSys_MATIS();

private:

    boost::shared_ptr<LocalToGlobalMap>   lg_map;
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
       PetscBool             pure_neumann;	/* PetscBool instead of PetscTruth*/
    } MatMyIS ;
};



#endif /* LA_LINSYS_HH_ */
