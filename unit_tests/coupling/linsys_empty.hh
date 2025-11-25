
#ifndef LINSYS_EMPTY_HH_
#define LA_LINSYS_PETSC_HH_

#include <functional>    // for unary_function
#include <string>        // for string
#include <vector>        // for vector
#include "la/linsys.hh"  // for LinSys
#include "petscksp.h"    // for KSP, KSPConvergedReason, _p_KSP
#include "petscmat.h"    // for Mat, MatCopy, MatZeroEntries, MatAssemblyType
#include "petscmath.h"   // for PetscScalar
#include "petscsys.h"    // for PetscErrorCode, PETSC_NULL
#include "petscvec.h"    // for Vec, _p_Vec, VecCopy, VecSet

class Distribution;

class LinSysEmpty : public LinSys
{

public:
    LinSysEmpty(const  Distribution * rows_ds)
    : LinSys( rows_ds ),
      matrix_(0),
	  rhs_(0) {}

    /**
     * Copy constructor.
     */
    LinSysEmpty( LinSys_PETSC &other )
    : LinSys(other),
      matrix_(0),
	  rhs_(0) {}


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

    PetscErrorCode mat_zero_entries() override {
        return 0;
    }

    PetscErrorCode rhs_zero_entries() override {
        return 0;
    }

    void start_allocation() override
    {}

    void start_add_assembly() override
    {}

    void mat_set_values( FMT_UNUSED int nrow, FMT_UNUSED int *rows, FMT_UNUSED int ncol, FMT_UNUSED int *cols, FMT_UNUSED double *vals ) override
    {}

    void rhs_set_values( FMT_UNUSED int nrow, FMT_UNUSED int *rows, FMT_UNUSED double *vals ) override
    {}

    void finish_assembly() override
    {}

    void set_tolerances(FMT_UNUSED double r_tol, FMT_UNUSED double a_tol, FMT_UNUSED double d_tol, FMT_UNUSED unsigned int max_it) override
    {}

    void apply_constrains( double scalar = 1. ) override
    {}

    LinSys::SolveInfo solve() override {
        return LinSys::SolveInfo(0, 0);
    }

    double get_solution_precision() override {
        return 0.0;
    }

    double compute_residual() override {
        return 0.0;
    }

    virtual ~LinSysEmpty()
    {}


protected:
    Mat     matrix_;             //!< Petsc matrix of the problem.
    Vec     rhs_;                //!< PETSc vector constructed with vx array.
};

#endif /* LINSYS_EMPTY_HH_ */
