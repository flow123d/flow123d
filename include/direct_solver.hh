/* 
 * File:   direct_solver.hh
 * Author: jb
 *
 * Created on May 19, 2010, 12:14 PM
 */

#ifndef _DIRECT_SOLVER_HH
#define	_DIRECT_SOLVER_HH

#include <lac/petsc_precondition.h>
#include <lac/petsc_solver.h>
#include <lac/petsc_sparse_matrix.h>
//#include <lac/petsc_block_sparse_matrix.h>

 #include <lac/petsc_vector.h>
 //#include <lac/petsc_parallel_vector.h>
 //#include <lac/petsc_parallel_sparse_matrix.h>

// TODO:
// DEAL sama poskytuje jen iteracni solvery
// rozhrani k PETSC, TRILONOs a UMFPacku, je nuten pouzit vcetne specialniho typu matic
// takze by bylo potreba mit sestavovani do obecne matice, aby to slo snadno menit


/**
 *  Complete solve by Schur complement
 */



class DirectSolver : public Subscriptor{



public:
    DirectSolver(PETScWrappers::SparseMatrix &m, PETScWrappers::Vector &rhs, PETScWrappers::Vector &x)
    :   system_matrix(m),
        system_rhs(rhs),
        system_x(x)
    {}

    void solve();

private:
    const PETScWrappers::SparseMatrix &system_matrix;
    const PETScWrappers::Vector &system_rhs;
    PETScWrappers::Vector &system_x;

};


void DirectSolver :: solve ()
{
   SolverControl solver_control(100,1.e-10);
   PETScWrappers :: SolverPreOnly solver(solver_control);
   PETScWrappers :: PreconditionLU precondLU(system_matrix);
   //PETScWrappers :: PreconditionerBase precondLU(system_matrix);
   //PETScWrappers :: MatrixBase mm;
   PETScWrappers :: Vector vv;
   solver.solve(system_matrix,system_x,system_rhs,precondLU);
   //solver.solve(system_matrix,vv,vv,precondLU);
}


#endif	/* _DIRECT_SOLVER_HH */

