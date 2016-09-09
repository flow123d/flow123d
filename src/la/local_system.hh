
#ifndef LOCAL_SYSTEM_HH_
#define LOCAL_SYSTEM_HH_

#include <armadillo>


/** Local system class is meant to be used for local assembly and then pass to global linear system.
 * The key idea is to take care of known solution values (Dirichlet boundary conditions) in a common way.
 * 
 * Usage the class consists of 3 steps:
 * 1) create local system, set global DoFs.
 * 2) set all known values (Dirichlet BC)
 * 3) set matrix and RHS entries 
 *    (if the entry is on dirichlet row or column, it is now taken care of)
 * 4) possibly fix the diagonal entries of the local system, where Dirichlet BC is set
 * 
 */
class LocalSystem
{
protected:
    arma::mat matrix;   ///< local system matrix
    arma::vec rhs;      ///< local system RHS
    
    /// vector of row indices where the solution is set (dirichlet BC)
    std::vector<unsigned int> loc_solution_rows;
    /// vector of solution values at @p loc_solution_rows indices (dirichlet BC)
    std::vector<double> loc_solution;
    
    std::vector<int> row_dofs;  ///< global row indices
    std::vector<int> col_dofs;  ///< global column indices
    
public:
    
    /** @brief Constructor.
     * 
     * @p nrows is number of rows of local system
     * @p ncols is number of columns of local system
     */
    LocalSystem(unsigned int nrows, unsigned int ncols);
    
    /// Resets the matrix, RHS, dofs to zero and clears solution settings
    void reset();
    
    /** @brief Set the position and value of known solution.
     * 
     * @p loc_row is local row index in solution vector
     * @p solution is the values of the solutioin
     */
    void set_solution(unsigned int loc_row, double solution);
    
    /** When the local system is assembled,
     * the diagonal entries on rows where the solution is set might be zero.
     * Therefore it is necessary to set a proper value to the diagonal entry
     * and respective RHS entry, such that the given solution holds.
     * 
     * Call only when the assembly of local system is finished.
     */
    void fix_diagonal();

    /**
    * This is a copy of the set_values method from LinSys.
    * 
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
    void set_value(unsigned int row_dof, unsigned int col_dof,
                    double mat_val, double rhs_val);        
    
    /**
        * This is a copy of the set_values method from LinSys.
        * 
    * Shortcut to assembly into matrix and RHS in one call, possibly apply Dirichlet boundary conditions.
    * @p row_dofs are global indices of rows of dense @p matrix and rows of dense vector @rhs in global system
    * @p col_dofs are global indices of columns of the matrix
    * @p loc_matrix is local matrix which is to be set in the local system
    * @p loc_rhs is local rhs which is to be set in the local system
    * @p row_solution are values of dofs belonging to row_dofs
    * @p col_solution are values of dofs belonging to col_dofs
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
    *   for passing local matrices into LinSys.
    *
    */
    void set_values(std::vector<int> &row_dofs, std::vector<int> &col_dofs,
                    const arma::mat &loc_matrix, const arma::vec &loc_rhs,
                    const arma::vec &row_solution, const arma::vec &col_solution);
};
    
#endif // LOCAL_SYSTEM_HH_