
#ifndef LOCAL_SYSTEM_HH_
#define LOCAL_SYSTEM_HH_

#include <armadillo>
#include "system/index_types.hh"

class LinSys;

/** Local system class is meant to be used for local assembly and then pass to global linear system.
 * The key idea is to take care of known solution values (Dirichlet boundary conditions) in a common way.
 * 
 * A] Usage of the class consists of 5 steps:
 * 1) create local system, set local/global DoFs.
 * 2) possibly setup sparsity pattern for the local system
 * 3) set matrix and RHS entries 
 * 4) set all known values (Dirichlet BC)
 *    (the order of items 3 and 4 does not matter)
 * 5) possibly eliminate the known solution
 *    (possibly fix the diagonal entries of the local system, where Dirichlet BC is set)
 * 6) pass the local system to the global system, possibly with the map local_to_global
 * 
 * B] Scenario with computation of Schur complement:
 * 1) do A1, do not have to set dof indices (if not doing B3)
 * 2) do A2, A3
 * 3) possibly do A4, A5, if there is known solution for the dofs not included in the Schur complement
 * 3) create another LocalSystem for Schur complement (local_schur), possibly set dof indices
 * 4) possibly do A2 on local_schur
 * 4) compute Schur complement, passing the prepared local_schur
 * 5) if not done before, set dof indices
 * 6) possibly do A4, A5 on local_schur
 * 7) do A6 with local_schur
 * 
 * C] Reconstruction from the Schur complement
 * 1) do A1, A3
 * 2) reconstruct the full solution
 * 
 * TODO:
 *  - set dofs as references
 *  - done: use arma::vec<int> to keep dofs (efficient, can be constructed from raw arrays)
 *  - done: use local dof indeces to set solution
 *  - rename set_solution to set_constraint
 *  -
 */
class LocalSystem
{
public:
    /**
     * Global row and col indices.  Are public and can be freely set.
     * Nevertheless one can also provide reference to already existing arrays through
     * specific constructor or reset function.
     */
    LocDofVec row_dofs, col_dofs;
    
    /**
     * @brief Default constructor.
     *
     * Object must be initialized by subsequent call of reset(nrows, ncols).
     */
    LocalSystem();

    /** @brief Constructor.
     * 
     * @p nrows is number of rows of local system
     * @p ncols is number of columns of local system
     */
    LocalSystem(uint nrows, uint ncols);
    
    /// Resets the matrix, RHS, dofs to zero and clears solution settings
    void reset();
    
    /// Resize and reset.
    void reset(uint nrows, uint ncols);

    /**
     * Resize and reset. Set dofs vectors to reuse arrays provided by given vectors.
     * Given vectors can not be changed until next call to of any reset function.
     */
    void reset(const LocDofVec &row_dofs, const LocDofVec &col_dofs);

    const arma::mat& get_matrix() {return matrix;}
    const arma::vec& get_rhs() {return rhs;}
    
    /** @brief Set the position and value of known solution. E.g. Dirichlet boundary condition.
     * 
     * @p loc_dofs is local row index in solution vector
     * @p solution is the values of the solution
     * @p diag_val is preferred diagonal value on the solution row
     */
    void set_solution(uint loc_dof, double solution, double diag=0.0);

    void set_solution_row(uint loc_row, double solution, double diag=0.0);

    void set_solution_col(uint loc_col, double solution);


    /** 
     * When finished with assembly of the local system,
     * this function eliminates all the known dofs.
     * 
     * It is skipped if there is not any solution dof set.
     * 
     * During elimination, the (global) diagonal entries on the rows, where the solution is set, might be zero.
     * Therefore it is necessary to set a proper value to the diagonal entry
     * and respective RHS entry, such that the given solution holds.
     * If preferred diagonal value has been set by @p set_solution then it is used.
     * 
     * Calling this function after the assembly of local system
     * and before passing the local system to the global one
     * is finished is users's responsibility.
     */
    void eliminate_solution();

    /** @brief Adds a single entry into the local system.
     * 
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     * @p rhs_val is RHS entry value
     */
    void add_value(uint row, uint col,
                   double mat_val, double rhs_val);
    
    /** @brief Matrix entry. 
     * Adds a single entry into the local system matrix.
     * 
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     */
    void add_value(uint row, uint col,
                   double mat_val);
    
    /** @brief RHS entry.
     * Adds a single entry into the local system RHS.
     * 
     * @p row is local row index of local system
     * @p rhs_val is RHS entry value
     */
    void add_value(uint row, double rhs_val);

    void set_matrix(arma::mat matrix);
    void set_rhs(arma::vec rhs);

    /// Sets the sparsity pattern for the local system.
    /** Due to petsc options: MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE)
     * all zeros will be thrown away from the system.
     * If we do not want some zero entries in the system matrix to be thrown away,
     * we can set these entries with this almost zero value.
     * 
     * Almost_zero values will be set in all entries: sp(i,j) != 0
     */
    void set_sparsity(const arma::umat & sp);
    
    /** @brief Computes Schur complement of the local system: S = C - B * invA * Bt
     * Applicable for square matrices.
     * It can be called either after eliminating Dirichlet dofs,
     * or the Dirichlet dofs can be set on the Schur complement
     * and the elimination done on the Schur complement.
     * 
     * @p offset index of the first row/column of submatrix C (size of A)
     * @p schur (output) LocalSystem with Schur complement
     * @p negative if true, the schur complement (including its rhs) is multiplied by -1.0
     */
    void compute_schur_complement(uint offset, LocalSystem& schur, bool negative=false) const;
    
    /** @brief Reconstructs the solution from the Schur complement solution: x = invA*b - invA * Bt * schur_solution
     * Applicable for square matrices.
     * 
     * @p offset index of the first row/column of submatrix C (size of A)
     * @p schur_solution solution of the Schur complement
     * @p reconstructed_solution (output) reconstructed solution of the complementary variable
     */
    void reconstruct_solution_schur(uint offset, const arma::vec &schur_solution, arma::vec& reconstructed_solution) const;

protected:
    void set_size(uint nrows, uint ncols);

    arma::mat matrix;   ///< local system matrix
    arma::vec rhs;      ///< local system RHS

    arma::mat sparsity; ///< sparsity pattern

    /// Number of rows/cols to be eliminated due to known solution.
    uint n_elim_rows, n_elim_cols;
    LocDofVec elim_rows;        /// Rows indices of rows to be eliminated.
    LocDofVec elim_cols;        /// Cols indices of cols to be eliminated.
    arma::vec solution_rows;    /// Values of the known solution (for row dofs).
    arma::vec solution_cols;    /// Values of the known solution (for col dofs).
    arma::vec diag_rows;        /// Prefered values on the diagonal after elimination.

    /// Due to petsc options: MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE)
    /// all zeros will be thrown away from the system.
    /// If we do not want some zero entries in the system matrix to be thrown away,
    /// we can set these entries with this almost zero value.
    ///
    /// This is done for example when BC values are eliminated and later the BC changes to different type (e.g. seepage).
    /// Another case is keeping the structure of matrix unchanged for the schur complements -
    /// for that we fill the whole diagonal (escpecially block C in darcy flow) with artificial zeros.
    static constexpr double almost_zero = std::numeric_limits<double>::min();


    friend class LinSys;
};
    
#endif // LOCAL_SYSTEM_HH_
