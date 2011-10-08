
#include <UMFPACK_system.hh>

/**
 * UMFPACK use only serial vectors
 * this is stright forward implementation
 *
 */ 

void UMFPACKSystem :: solve (VectorType &solution)
{
    condense_constraints();

    SparseDirectUMFPACK solver;
    Vector<double> RHS;
    RHS = system_rhs;
    Vector<double> sol(RHS);


   /*system_matrix.block(0,0).print(cout);
   system_matrix.block(0,1).print(cout);
   system_matrix.block(1,0).print(cout);
   system_matrix.block(1,1).print(cout);
   RHS.print(cout);
*/
    solver.initialize(system_matrix);
    solver.vmult(sol,RHS);
    solution=sol;
    
    constraints.distribute(solution);
    /*
    unsigned int i_sol=0;
    for(unsigned int i=0;i<active_blocks;i++) 
        for(unsigned int j=0;j<solution.block(i).size();j++) solution.block(i)(j) = sol(i_sol);
     */
}

void UMFPACKSystem::eliminate_block(VectorType &sol,const unsigned int block) 
{
    condense_constraints();

    SparseDirectUMFPACK solver;
    Vector<double> rhs(sol.block(block));

    // substitute pressure into velocity equation
    rhs=0;
    for(unsigned int i=0;i<blocks.size();i++)
        if (i != block)
            system_matrix.block(block,i).vmult_add (rhs, sol.block(i));
    rhs *= -1;
    rhs+=system_rhs.block(block);

    solver.initialize(system_matrix.block(block,block));
    solver.vmult(sol.block(block),rhs);
    constraints.distribute(sol);
}

//! compute residuum for given solution vector
double UMFPACKSystem::residual(const VectorType &sol, VectorType &res)
{

    condense_constraints();
    /* this was my own WRONG implementation of block residual
     * it was necessary for solving on partial block vectors
     *
    double i_norm;
    double norm=0;
    for(unsigned int i=0;i<blocks.size();i++)
    {
        res.block(i)=system_rhs.block(i);
        for(unsigned int j=0; j<blocks.size();j++) {
            system_matrix.block(i,j).vmult_add(res.block(i),sol.block(j));
        }
        i_norm=res.block(i).l2_norm();
        cout<<"block norm "<<i<<" = "<<i_norm<<endl;
        norm+=i_norm*i_norm;
    }
     */
    return system_matrix.residual(res,sol,system_rhs);
    
 }