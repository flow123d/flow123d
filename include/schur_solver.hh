/*
 * File:   schur_solver.hh
 * Author: jb
 *
 * Created on May 17, 2010, 5:43 PM
 */

#ifndef _SCHUR_SOLVER_HH
#define	_SCHUR_SOLVER_HH

#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/block_vector.h>
#include <lac/sparse_direct.h>
#include <lac/sparse_ilu.h>
using namespace dealii;

#include <boost/assign.hpp>
using namespace boost::assign;

/**
 *  Implicit inverse for a given matrix.
 *  - use an ini vector as an solution approximation
 */

template <class MAT, class PRE = PreconditionIdentity>
class ImplicitInverse : public Subscriptor
{
  public:
    ImplicitInverse (const MAT &A, Vector<double> &ini_vector, const PRE &pre = PreconditionIdentity(), double tol = 1.0e-2)
            : matrix(&A), ini_vector(&ini_vector), r_tol(tol), precond(&pre)
    {}

    void vmult (Vector<double>  &dst,const Vector<double> &src) const
    {
        //cout << "------ internal solver " << r_tol <<endl;
        SolverControl solver_control (10000, r_tol * src.l2_norm());
        SolverCG<>    cg (solver_control);
        //if (&dst != ini_vector) dst = *ini_vector;
        //dst=0;
        cg.solve (*matrix, dst, src, *precond);
        //if (&dst != ini_vector) *ini_vector = dst;
    }

  private:
    const SmartPointer<const MAT > matrix;
    const SmartPointer<Vector<double> > ini_vector;
    double r_tol;
    const SmartPointer<const PRE > precond;
    
};


/**
 *  Full implicit Schur complement matrix
 *  - solution of inverse approximated by last solution
 */


class SchurComplement : public Subscriptor
{
    public:
        SchurComplement (const BlockSparseMatrix<double> &A, SparseILU<double> &precond)
                :
                IA_B_x_b (A.block(0,0).m()),
                system_matrix (&A),
                inner_preconditioner(precond),
                B_x_b (A.block(0,0).m()),
                C_x_b (A.block(1,1).m())
        {
            IA_B_x_b=0;

        }

        void vmult (Vector<double> &dst,  const Vector<double> &src) const
        {
            system_matrix->block(0,1).vmult (B_x_b, src);  // B * x_b

            ImplicitInverse<SparseMatrix<double>, SparseILU<double> >
                IA(system_matrix->block(0,0), IA_B_x_b, inner_preconditioner, 1e-6 );
            IA.vmult(IA_B_x_b, B_x_b);

            system_matrix->block(1,0).vmult (dst, IA_B_x_b);
            //dst*=-1;
            system_matrix->block(1,1).vmult(C_x_b,src);
            dst-=C_x_b;
        }
        mutable Vector<double> IA_B_x_b;

    private:
        const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
        const SparseILU<double> &inner_preconditioner;
        mutable Vector<double> B_x_b, C_x_b;
};


class SchurComplementApprox : public Subscriptor
{
    public:
        SchurComplementApprox (const BlockSparseMatrix<double> &A)
                :
                IA_B_x_b (A.block(0,0).m()),
                B_x_b (A.block(0,0).m()),
                C_x_b (A.block(1,1).m()),
                system_matrix (&A)
        {
            IA_B_x_b=0;

        }

        void vmult (Vector<double> &dst,  const Vector<double> &src) const
        {
            system_matrix->block(0,1).vmult (B_x_b, src);  // B * x_b
            system_matrix->block(0,0).precondition_Jacobi (IA_B_x_b, B_x_b);
            system_matrix->block(1,0).vmult (dst, IA_B_x_b);
            //dst*=-1;
            system_matrix->block(1,1).vmult(C_x_b,src);
            dst-=C_x_b;
        }
        mutable Vector<double> IA_B_x_b;

    private:
        mutable Vector<double> B_x_b, C_x_b;
        const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
};



/**
 *  Jacobi Precond matrix
 */

class JacobiPrecondMatrix : public Subscriptor
{
  public:
    JacobiPrecondMatrix(const SparseMatrix<double> &A, Vector<double> &v,  PreconditionIdentity p, double tol) : matrix(&A) {};

    void vmult (Vector<double>  &dst, const Vector<double> &src) const
    {
        matrix->precondition_Jacobi (dst, src);
    }

  private:
    const SmartPointer<const SparseMatrix<double> > matrix;
};

/** ***********************************
 *  Complete solve by Schur complement
 */
template <class DF>
class SchurSystemImplicit : public Subscriptor
{
public:
    typedef  BlockVector<double>        VectorType;
    typedef  BlockSparseMatrix<double>  MatrixType;

    SchurSystemImplicit(DF &dof_handler);
    void solve();
    inline void matrix_add(const std::vector< unsigned int > &indices, const FullMatrix< double > &full_matrix)
        { system_matrix.add(indices,full_matrix); }
    inline void rhs_add(const std::vector< unsigned int > &indices, const Vector<double>  &full_vec)
        { for (unsigned int i=0; i<indices.size(); ++i) system_rhs(indices[i]) += full_vec(i); }

    inline MatrixType &get_matrix() { return system_matrix; }
    inline VectorType &get_rhs() { return system_rhs; }
    inline VectorType &get_solution() { return solution; }

private:
    BlockSparsityPattern sparsity_pattern;
    SparseILU<double> inner_preconditioner;

    MatrixType system_matrix;
    VectorType system_rhs;
    VectorType solution;

    int size_a;     // size of eliminated block A
    int size_b;     // size of the Schur complement


    Vector<double> IA_f_a;
    Vector<double> tmp_a, tmp_b, tmp_pre;
    Vector<double> schur_rhs;

    double r_tol;

};

template <class DF>
SchurSystemImplicit<DF> :: SchurSystemImplicit(DF &dof_handler)
{
    int dim=dof_handler.dimension;
    // renuber DoFs to get block global matrix
    std::vector<unsigned int> target_component (dim+1,1);
    target_component[dim] = 0;
  
   DoFRenumbering::component_wise (dof_handler,target_component);

   // get DOF handler parameters

   std::vector<unsigned int> dofs_per_component (dim+1);
   DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
   size_a = dofs_per_component[dim],
   size_b = dofs_per_component[0];
   std::cout<< "blocks: "<<size_a << size_b<< std::endl;
   const unsigned int
     n_couplings = dof_handler.max_couplings_between_dofs();

   // create global matrix

   sparsity_pattern.reinit (2,2);
   sparsity_pattern.block(0,0).reinit (size_a, size_a, n_couplings);
   sparsity_pattern.block(1,0).reinit (size_b, size_a, n_couplings);
   sparsity_pattern.block(0,1).reinit (size_a, size_b, n_couplings);
   sparsity_pattern.block(1,1).reinit (size_b, size_b, n_couplings);
   sparsity_pattern.collect_sizes();

   // sparsity pattern is given by interrelation of DoFHandler, FE, and Triangulation
   DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
   sparsity_pattern.compress();

   system_matrix.reinit (sparsity_pattern);

   // create Solution and RHS vectors
   solution.reinit (2);
   solution.block(0).reinit (size_a);
   solution.block(1).reinit (size_b);
   solution.collect_sizes ();

   system_rhs.reinit ( solution );

   // Create Schur solver

   IA_f_a.reinit(size_a);
   tmp_a.reinit(size_a);
   tmp_b.reinit(size_b);
   tmp_pre.reinit(size_b);
   schur_rhs.reinit(size_b);


    IA_f_a=0;
    tmp_pre=0;
    r_tol=1e-8;
}



template <class DF>
void SchurSystemImplicit<DF> :: solve ()
{
    inner_preconditioner.initialize (system_matrix.block(0,0),
       SparseILU<double>::AdditionalData());


  //std::ofstream out ("sparsity_pattern.gpl");
  //sparsity_pattern.print_gnuplot(out);

  //system_matrix.compress();

  {
    // form Schur RHS:  Bt IA f_a - f_b
    // need an approximation of IA f_a
    // use the last solution
    
    
    ImplicitInverse<SparseMatrix<double>, SparseILU<double> >
        IA(system_matrix.block(0,0), IA_f_a,inner_preconditioner, r_tol);
    IA.vmult(IA_f_a, system_rhs.block(0));

    system_matrix.block(1,0).vmult (schur_rhs, IA_f_a);
    schur_rhs -= system_rhs.block(1);
    //schur_rhs -= tmp_b;
    //schur_rhs *=-1;
  }

    SchurComplement schur_complement(system_matrix, inner_preconditioner);
    SchurComplementApprox schur_complement_approx(system_matrix);
    //SchurComplement<JacobiPrecondMatrix>                    schur_complement_jacobi(system_matrix);
    //PreconditionJacobi<SparseMatrix<double> >               precond_jacobi;
    ImplicitInverse<SchurComplementApprox> schur_complement_precond(schur_complement_approx, tmp_pre, PreconditionIdentity(),  1e-2);
    //, PreconditionJacobi<SparseMatrix<double> > >
    //    schur_complement_precond(schur_complement_jacobi, tmp_pre, precond_jacobi,  1e-2);

    //precond_jacobi.initialize( system_matrix.block(1,1) );


  {
  // Main CG solver
  SolverControl solver_control ( size_b, r_tol * schur_rhs.l2_norm());
  solver_control.log_history(true);
  solver_control.log_frequency(1);
  SolverCG<>    cg (solver_control);

  //PreconditionJacobi<SparseMatrix<double> > precond;
  //precond.initialize( system_matrix.block(1,1) );
  cg.solve (schur_complement, solution.block(1), schur_rhs,
            //precond) ;
            schur_complement_precond);
      std::cout << "iter: "<< solver_control.last_step()<<std::endl;
  }

  // recover x_a = IA * f_a - IA * B * x_b
  {
    system_matrix.block(0,1).vmult (tmp_a, solution.block(1));
    tmp_a *= -1;
    tmp_a += system_rhs.block(0);
    solution.block(0) = IA_f_a;
    solution.block(0) -= schur_complement.IA_B_x_b;     // solution approx.

    ImplicitInverse<SparseMatrix<double>, SparseILU<double> >
        IA(system_matrix.block(0,0), solution.block(0),inner_preconditioner, r_tol);
    IA.vmult(solution.block(0), tmp_a);
  }
}
#endif	/* _SCHUR_SOLVER_HH */

