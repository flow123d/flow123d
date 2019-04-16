
#ifndef LOCAL_SYSTEM_HH_
#define LOCAL_SYSTEM_HH_

#include <armadillo>
#include "mesh/long_idx.hh"

class LinSys;

/** Local system class is meant to be used for local assembly and then pass to global linear system.
 * The key idea is to take care of known solution values (Dirichlet boundary conditions) in a common way.
 * 
 * Usage of the class consists of 3 steps:
 * 1) create local system, set global DoFs.
 * 2) set all known values (Dirichlet BC)
 * 3) set matrix and RHS entries 
 *    (if the entry is on dirichlet row or column, it is now taken care of)
 * 4) eliminate known solution and possibly fix the diagonal entries of the local system, where Dirichlet BC is set
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
    typedef arma::Col<LongIdx> DofVec;

    /**
     * Global row and col indices.  Are public and can be freely set.
     * Nevertheless one can also provide reference to already existing arrays through
     * specific constructor or reset function.
     */
    DofVec row_dofs, col_dofs;
    
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
    LocalSystem(unsigned int nrows, unsigned int ncols);
    
    /// Resets the matrix, RHS, dofs to zero and clears solution settings
    void reset();
    
    /// Resize and reset.
    void reset(arma::uword nrows, arma::uword ncols);

    /**
     * Resize and reset. Set dofs vectors to reuse arrays provided by given vectors.
     * Given vectors can not be changed until next call to of any reset function.
     */
    void reset(const DofVec &row_dofs, const DofVec &col_dofs);

    const arma::mat& get_matrix() {return matrix;}
    const arma::vec& get_rhs() {return rhs;}
    
    /** @brief Set the position and value of known solution. E.g. Dirichlet boundary condition.
     * 
     * @p loc_dofs is local row index in solution vector
     * @p solution is the values of the solution
     * @p diag_val is preferred diagonal value on the solution row
     */
    void set_solution(unsigned int loc_dof, double solution, double diag=0.0);

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
     * Calling this function after the assembly of local system is finished is users's responsibility.
     */
    void eliminate_solution();

    /** @brief Adds a single entry into the local system.
     * 
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     * @p rhs_val is RHS entry value
     */
    void add_value(unsigned int row, unsigned int col,
                   double mat_val, double rhs_val);
    
    /** @brief Matrix entry. 
     * Adds a single entry into the local system matrix.
     * 
     * @p row is local row index of local system
     * @p col is local column index of local system
     * @p mat_val is matrix entry value
     */
    void add_value(unsigned int row, unsigned int col,
                   double mat_val);
    
    /** @brief RHS entry.
     * Adds a single entry into the local system RHS.
     * 
     * @p row is local row index of local system
     * @p rhs_val is RHS entry value
     */
    void add_value(unsigned int row, double rhs_val);

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
     * 
     * @p offset index of the first row/column of submatrix C (size of A)
     * @p schur (output) LocalSystem with Schur complement
     */
    void compute_schur_complement(uint offset, LocalSystem& schur);
    
    /** @brief Reconstructs the solution from the Schur complement solution: x = invA*b - invA * Bt * schur_solution
     * Applicable for square matrices.
     * 
     * @p offset index of the first row/column of submatrix C (size of A)
     * @p schur_solution solution of the Schur complement
     * @p reconstructed_solution (output) reconstructed solution of the complementary variable
     */
    void reconstruct_solution_schur(uint offset, const arma::vec &schur_solution, arma::vec& reconstructed_solution);

protected:
    void set_size(unsigned int nrows, unsigned int ncols);

    arma::mat matrix;   ///< local system matrix
    arma::vec rhs;      ///< local system RHS

    arma::mat sparsity; ///< sparsity pattern

    unsigned int n_elim_rows, n_elim_cols;
    DofVec elim_rows;
    DofVec elim_cols;
    arma::vec solution_rows;
    arma::vec solution_cols;
    arma::vec diag_rows;

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
