/* 
 * File:   UMFPACK_solver.hh
 * Author: jb
 *
 * Created on May 26, 2010, 6:36 PM
 */

#ifndef _UMFPACK_SOLVER_HH
#define	_UMFPACK_SOLVER_HH

#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/block_vector.h>
#include <lac/sparse_direct.h>


#include <boost/assign.hpp>
using namespace boost::assign;

/** ***********************************
 *  Direct solve of block matrix
 */
template <class DF>
class UMFPACKSystem : public Subscriptor
{
public:
    typedef  BlockVector<double>        VectorType;
    typedef  BlockSparseMatrix<double>  MatrixType;

    UMFPACKSystem(DF &dof_handler);
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

    MatrixType system_matrix;
    VectorType system_rhs;
    VectorType solution;
    
    

    int size_a;     // size of eliminated block A
    int size_b;     // size of the Schur complement

};

template <class DF>
UMFPACKSystem<DF> :: UMFPACKSystem(DF &dof_handler)
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
   cout<< "blocks: "<<size_a << size_b<< endl;
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

}



template <class DF>
void UMFPACKSystem<DF> :: solve ()
{
    SparseDirectUMFPACK solver;
    Vector<double> RHS;
    RHS = system_rhs;
    Vector<double> sol(RHS);

    solver.initialize(system_matrix);
    solver.vmult(sol,RHS);
    solution = sol;

    BlockVector<double> resid(solution);
    system_matrix.vmult(resid,solution);
    resid-=system_rhs;

    //cout<< "Residuum: "<< resid.l2_norm()<<endl;
}

#endif	/* _UMFPACK_SOLVER_HH */

