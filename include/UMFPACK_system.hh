/* 
 * File:   UMFPACK_solver.hh
 * Author: jb
 *
 * Created on May 26, 2010, 6:36 PM
 */
/**
 * Aim is to have an abstract class for block linear problem, that
 * provides global assembly operations, method for resolving velocities
 * and block access to vectors. The solution vector is outside of the object
 * and can have more blocks then those used for solve.
 *
 * Problem:
 * Chci resit Richadse ruznymi solvery, k tomu potrebuju
 * 1) metody lin. seystemu pro sestavovani (matici a rhs me nezajima)
 * 2) metodu solv, ktera bude pracovat blokovym vektorem
 *    - bud bude solver vyzadovat nejaky typ blokoveho vektoru a
 *      ten bude dostupny venku a bude poskytovat potrebne oprace
 *      a to jak pro venek, tak pro vnitrek !!!
 *
 * pres templaty: nemam zadnou virtualni tridu pro jednotlive solvery, ale
 * umim zkontrolovat, ze implementuji dane metody a maji typ blokoveho vektoru,
 * ktery poskytuje pozadovane operace (?? co venku potrebuju delat v blokovym vektorem?
 * asi zejmena pristup k jednotlivym blokum, linearni pocitani s vektory a normy
 *
 * pres virtualni tridu to nejde
 * rozhrani virtualni tridy uz musi predepsat typ blokoveho vektoru
 *
 * varianta, trida templatovana typem vektoru a se specializacemi
 *
 * zatim mam jen UMFPACKimplementaci a predpokladam dalsi rozvoj podle varianty template
 *
 * BlockVectorBase je spolecnym predkem vsech Block vektoru, to by slo vyuzit,
 * dale by trida mela byt zamerena na reseni pomoci Schur doplnku.
 * Asi nejlepsi by bylo mit implmentaci s PETSC vektorama a mit soucasnou zjednodusujici variantu s UMFPAckem
 * 
 */


#ifndef _UMFPACK_SOLVER_HH
#define	_UMFPACK_SOLVER_HH

#include <fstream>
#include <iostream>


#include <base/subscriptor.h>
#include <lac/block_sparse_matrix.h>
#include <lac/precondition.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/block_vector.h>
#include <lac/sparse_direct.h>

#include <dofs/dof_tools.h>
using namespace dealii;
using namespace std;
//#include <boost/assign.hpp>
//using namespace boost::assign;

/** ***********************************
 *  Direct solve of block matrix
 */
class UMFPACKSystem : public Subscriptor
{
public:
    typedef  BlockVector<double>        VectorType;
    typedef  BlockSparseMatrix<double>  MatrixType;
    typedef  BlockSparsityPattern       PatternType;

    UMFPACKSystem(unsigned int n_blocks)
        : blocks(n_blocks) {}
    
    template <class DH>
    void reinit(const DH &dof_handler);
    //! solve linear system using the acitve blocks
    void solve(VectorType &sol);
    //! elimiate given block
    void eliminate_block(VectorType &sol,unsigned int block);
    //! compute residuum for given solution vector
    double residual(const VectorType &sol, VectorType &res);

    inline void condense_constraints()
    {
        if (! condensed) {
        constraints.condense(system_matrix);
        constraints.condense(system_rhs);
        condensed=true;
        }
    }

    //! assembly methods
    inline void set_zero()
        { system_matrix=0; system_rhs=0; condensed=false; }
    inline void matrix_add(const std::vector< unsigned int > &indices, const FullMatrix< double > &full_matrix)
        { system_matrix.add(indices,full_matrix); }
    inline void rhs_add(const std::vector< unsigned int > &indices, const Vector<double>  &full_vec)
        { for (unsigned int i=0; i<indices.size(); ++i) system_rhs(indices[i]) += full_vec(i); }

    //! exceptional access methods, these shloud be somehow eliminated
    inline MatrixType &get_matrix() { return system_matrix; }
    inline VectorType &get_rhs() { return system_rhs; }
    inline PatternType &get_pattern() { return sparsity_pattern; }
    inline ConstraintMatrix &get_constraints() { return constraints; }

private:
    std::vector<unsigned int> blocks;
    mutable bool condensed;

    PatternType sparsity_pattern;
    MatrixType system_matrix;
    VectorType system_rhs;

    ConstraintMatrix constraints;
};


template <class DH>
void UMFPACKSystem :: reinit(const DH &dof_handler)
{
   constraints.clear();
   DoFTools::make_hanging_node_constraints (dof_handler, constraints);
   constraints.close();
   DoFTools::count_dofs_per_block (dof_handler, blocks);

   
   std::cout<< "blocks: ";
   copy(blocks.begin(), blocks.end(), ostream_iterator<int>( cout, " " ) );
   cout << std::endl;

   const unsigned int n_couplings = dof_handler.max_couplings_between_dofs();

   // create global matrix
   system_matrix.clear();
   sparsity_pattern.reinit (blocks.size(),blocks.size());
   //system_matrix.reinit (n_blocks,n_blocks);
   for(unsigned int i=0;i<blocks.size();i++)
       for(unsigned int j=0;j<blocks.size();j++)
            sparsity_pattern.block(i,j).reinit (blocks[i], blocks[j], n_couplings);
   sparsity_pattern.collect_sizes();
   DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern, constraints);
   constraints.condense(sparsity_pattern);
   sparsity_pattern.compress();

   
   system_matrix.reinit (sparsity_pattern);

    system_rhs.reinit(blocks.size());
    for (unsigned int i=0; i < blocks.size(); i++) {
        system_rhs.block(i).reinit(blocks[i]); // fast reinit, since
    }
    system_rhs.collect_sizes();

    condensed=false;
        
}

#endif	/* _UMFPACK_SOLVER_HH */

