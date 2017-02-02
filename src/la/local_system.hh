
#ifndef LOCAL_SYSTEM_HH_
#define LOCAL_SYSTEM_HH_

#include <armadillo>

class LinSys;

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
friend class LinSys;
protected:
    arma::mat matrix;   ///< local system matrix
    arma::vec rhs;      ///< local system RHS
    
    /// vector of row indices where the solution is set (dirichlet BC)
    std::vector<unsigned int> loc_solution_rows;
    /// vector of solution values at @p loc_solution_rows indices (dirichlet BC)
    std::vector<double> loc_solution;
    /// diagonal values for dirichlet BC rows (set in set_value)
    std::vector<double> diag_values;
    /// diagonal values for dirichlet BC rows (set in set_solution)
    std::vector<double> preferred_diag_values;
    
public:
    
    std::vector<int> row_dofs;  ///< global row indices
    std::vector<int> col_dofs;  ///< global column indices
    
    /** @brief Constructor.
     * 
     * @p nrows is number of rows of local system
     * @p ncols is number of columns of local system
     */
    LocalSystem(unsigned int nrows, unsigned int ncols);
    
    /// Resets the matrix, RHS, dofs to zero and clears solution settings
    void reset();
    
    const arma::mat& get_matrix() {return matrix;}
    const arma::vec& get_rhs() {return rhs;}
    
    /** @brief Set the position and value of known solution.
     * 
     * @p loc_row is local row index in solution vector
     * @p solution is the values of the solutioin
     */
    void set_solution(unsigned int loc_row, double solution, double diag_val = 0.0);
    
    /** When the local system is assembled,
     * the diagonal entries on rows, where the solution is set, might be zero.
     * Therefore it is necessary to set a proper value to the diagonal entry
     * and respective RHS entry, such that the given solution holds.
     * 
     * Call only when the assembly of local system is finished.
     */
    void fix_diagonal();

    /** @brief Sets a single entry into the local system.
     * 
     * Known solution must be set before, so it is eliminated correctly during his call.
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     * @p rhs_val is RHS entry value
     */
    void set_value(unsigned int row, unsigned int col,
                   double mat_val, double rhs_val);

    /** @brief Adds a single entry into the local system.
     * 
     * Known solution must be set before, so it is eliminated correctly during his call.
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     * @p rhs_val is RHS entry value
     */
    void add_value(unsigned int row, unsigned int col,
                   double mat_val, double rhs_val);
    
//     /** @brief Sets a single entry into the local system matrix.
//      * 
//      * Known solution must be set before, so it is eliminated correctly during his call.
//      * @p row is local row index of local system
//      * @p col is local column index of local system
//      * @p mat_val is matrix entry value
//      */
//     void set_mat_value(unsigned int row, unsigned int col, double mat_val);
    
    /** @brief Sets a submatrix and rhs subvector into the local system matrix.
     * 
     * Known solution must be set before, so it is eliminated correctly during his call.
     * @p rows are local row indices of local system
     * @p cols are local column indices of local system
     * @p loc_mat is submatrix to be entered
     */
    void set_mat_values(const std::vector< unsigned int >& rows,
                        const std::vector< unsigned int >& cols,
                        const arma::mat& loc_matrix);
    
    /** @brief Sets a submatrix and rhs subvector into the local system.
     * 
     * Known solution must be set before, so it is eliminated correctly during his call.
     * @p rows are local row indices of local system
     * @p cols are local column indices of local system
     * @p loc_mat is submatrix to be entered
     * @p loc_rhs is vector to entered into RHS
     */
    void set_values(const std::vector<unsigned int> &rows,
                    const std::vector<unsigned int> &cols,
                    const arma::mat &loc_matrix,
                    const arma::vec &loc_rhs);
    
    /**
     * This is a copy of the set_values method from LinSys. (OBSOLETE)
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
//     void set_values(std::vector<int> &rows, std::vector<int> &cols,
//                     const arma::mat &loc_matrix, const arma::vec &loc_rhs,
//                     const arma::vec &row_solution, const arma::vec &col_solution);
};
    
#endif // LOCAL_SYSTEM_HH_