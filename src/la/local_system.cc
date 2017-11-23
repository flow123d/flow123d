
#include "local_system.hh"

#include <armadillo>
#include "system/sys_vector.hh"


LocalSystem::LocalSystem()
{}


LocalSystem::LocalSystem(unsigned int nrows, unsigned int ncols)
: row_dofs(nrows),
  col_dofs(ncols),
  matrix(nrows, ncols),
  rhs(nrows),
  sparsity_pattern_(nrows,ncols),
  elim_rows(nrows),
  elim_cols(ncols),
  solution_rows(nrows),
  solution_cols(ncols),
  diag_rows(nrows)
{
    reset();
}

void LocalSystem::set_size(unsigned int nrows, unsigned int ncols)
{
    row_dofs.set_size(nrows);
    col_dofs.set_size(ncols);
    matrix.set_size(nrows, ncols);
    rhs.set_size(nrows);
    elim_rows.set_size(nrows);
    elim_cols.set_size(ncols);
    solution_rows.set_size(nrows);
    solution_cols.set_size(ncols);
    diag_rows.set_size(nrows);
    sparsity_pattern_.resize(nrows,ncols);
}


void LocalSystem::reset()
{
    // zeros in local system
    matrix.zeros();
    rhs.zeros();
    // drop all dirichlet values
    n_elim_rows=n_elim_cols=0;
}


void LocalSystem::reset(unsigned int nrows, unsigned int ncols)
{
    matrix.set_size(nrows, ncols);
    rhs.set_size(nrows);
    row_dofs.resize(matrix.n_rows);
    col_dofs.resize(matrix.n_cols);
    sparsity_pattern_.resize(nrows,ncols);
    reset();
}



void LocalSystem::reset(const DofVec &rdofs, const DofVec &cdofs)
{
    set_size(rdofs.n_rows, cdofs.n_rows);
    sparsity_pattern_.resize(rdofs.n_rows,cdofs.n_rows);
    reset();
    row_dofs = rdofs;
    col_dofs = cdofs;
}



void LocalSystem::set_solution(unsigned int loc_dof, double solution, double diag)
{
    // check that dofs are same
    //ASSERT_DBG( arma::all(row_dofs == col_dofs) );
    set_solution_row(loc_dof, solution, diag);
    set_solution_col(loc_dof, solution);
}

void LocalSystem::set_solution_row(uint loc_row, double solution, double diag) {
    elim_rows[n_elim_rows]=loc_row;
    solution_rows[n_elim_rows] = solution;
    diag_rows[n_elim_rows] = diag;
    n_elim_rows++;
}

void LocalSystem::set_solution_col(uint loc_col, double solution) {
    elim_cols[n_elim_cols]=loc_col;
    solution_cols[n_elim_cols] = solution;
    n_elim_cols++;
}

/*

void LocalSystem::set_solution(const DofVec & loc_rows, const arma::vec &solution, const arma::vec &diag) {
    ASSERT_EQ_DBG(loc_rows.n_rows(), solution)
    set_solution_rows(loc_rows, solution, diag);
    set_solution_cols()
}
void LocalSystem::set_solution_rows(DofVec & loc_rows, const arma::vec &solution, const arma::vec &diag);
void LocalSystem::set_solution_cols(DofVec & loc_cols, const arma::vec &solution);

*/

void LocalSystem::eliminate_solution()
{
    if (! n_elim_rows &&  ! n_elim_cols ) return;
    
    //DebugOut().fmt("elim rows: {} elim_cols: {}", n_elim_rows, n_elim_cols);
    
    arma::mat tmp_mat = matrix;
    arma::vec tmp_rhs = rhs;
    
    unsigned int ic, ir, row, col, j;

    // eliminate columns
    for(ic=0; ic < n_elim_cols; ic++) {
        col = elim_cols[ic];
        tmp_rhs -= solution_cols[ic] * tmp_mat.col( col );
        
        // keep the structure of the matrix by setting almost_zero when eliminating dirichlet BC
        for(j=0; j < tmp_mat.n_rows; j++)
            if(tmp_mat.col(col)(j) != 0) tmp_mat.col(col)(j) = almost_zero;
    }

    // eliminate rows
    for(ir=0; ir < n_elim_rows; ir++) {
        row = elim_rows[ir];
        tmp_rhs( row ) = 0.0;
        
        // keep the structure of the matrix by setting almost_zero when eliminating dirichlet BC
        for(j=0; j < tmp_mat.n_cols; j++)
            if(tmp_mat.row(row)(j) != 0) tmp_mat.row(row)(j) = almost_zero;

        // fix global diagonal
        for(ic=0; ic < n_elim_cols; ic++) {
            col = elim_cols[ic];
            if (row_dofs[row] == col_dofs[col]) {
                ASSERT_DBG(fabs(solution_rows[ir] - solution_cols[ic]) <1e-12 );
                // if preferred value is not set, then try using matrix value
                double new_diagonal = matrix(row, col);

                if (diag_rows[ir] != 0.0)    // if preferred value is set
                    new_diagonal = diag_rows[ir];
                else if(new_diagonal == 0.0)      // if an assembled value is not available
                    new_diagonal = 1.0;
                //                     double new_diagonal = fabs(matrix(sol_row,col));
                //                     if (new_diagonal == 0.0) {
                //                         if (matrix.is_square()) {
                //                             new_diagonal = arma::sum( abs(matrix.diag())) / matrix.n_rows;
                //                         } else {
                //                             new_diagonal = arma::accu( abs(matrix) ) / matrix.n_elem;
                //                         }
                //                     }
                tmp_mat(row,col) = new_diagonal;
                tmp_rhs(row) = new_diagonal * solution_rows[ir];

            }
        }
    }

    matrix = tmp_mat;
    
    // filling almost_zero according to sparsity pattern
    for(unsigned int i=0; i < matrix.n_rows; i++)
        for(unsigned int j=0; j < matrix.n_cols; j++)
        {
            // Due to petsc options: MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE)
            // all zeros will be thrown away from the system.
            // Since we need to keep the matrix structure for the schur complements
            // and time dependent flow, we fill the whole diagonal with artificial zeros.
            if (sparsity_pattern_.is_set_global_diagonal() && row_dofs[i] == col_dofs[j])
                matrix(i,j) = matrix(i,j) + almost_zero;
            else
                matrix(i,j) = matrix(i,j) + almost_zero * sparsity_pattern_.get_matrix()(i,j);
        }
        
    rhs = tmp_rhs;
    n_elim_cols=n_elim_rows=0;

    //DebugOut() << matrix;
    //DebugOut() << rhs;


}


void LocalSystem::add_value(unsigned int row, unsigned int col, double mat_val, double rhs_val)
{
    ASSERT_DBG(row < matrix.n_rows);
    ASSERT_DBG(col < matrix.n_cols);
    
    matrix(row, col) += mat_val;
    rhs(row) += rhs_val;
}

void LocalSystem::add_value(unsigned int row, unsigned int col, double mat_val)
{
    ASSERT_DBG(row < matrix.n_rows);
    ASSERT_DBG(col < matrix.n_cols);
    
    matrix(row, col) += mat_val;
}

void LocalSystem::add_value(unsigned int row, double rhs_val)
{
    ASSERT_DBG(row < matrix.n_rows);
    
    rhs(row) += rhs_val;
}


void LocalSystem::set_matrix(arma::mat m) {
    ASSERT_EQ_DBG(matrix.n_rows, m.n_rows);
    ASSERT_EQ_DBG(matrix.n_cols, m.n_cols);
    matrix = m;
}

void LocalSystem::set_rhs(arma::vec r) {
    ASSERT_EQ_DBG(matrix.n_rows, r.n_rows);
    rhs = r;
}

LocalSystem::SparsityPattern & LocalSystem::sparsity_pattern()
{
    return sparsity_pattern_;
}

LocalSystem::SparsityPattern::SparsityPattern()
  : force_global_diagonal_(false)
{}

LocalSystem::SparsityPattern::SparsityPattern(unsigned int nrows, unsigned int ncols)
  : force_global_diagonal_(false)
{
    pattern.zeros(nrows,ncols);
}

void LocalSystem::SparsityPattern::set(unsigned int row, unsigned int col)
{
    ASSERT_DBG(row < pattern.n_rows);
    ASSERT_DBG(col < pattern.n_cols);
    pattern(row,col) = 1;
}

void LocalSystem::SparsityPattern::reset(unsigned int row, unsigned int col)
{
    ASSERT_DBG(row < pattern.n_rows);
    ASSERT_DBG(col < pattern.n_cols);
    pattern(row,col) = 0;
}

void LocalSystem::SparsityPattern::fill_all()
{
    pattern.ones();
}

void LocalSystem::SparsityPattern::clear_all()
{
    pattern.zeros();
}


void LocalSystem::SparsityPattern::force_global_diagonal()
{
    force_global_diagonal_ = true;
}

bool LocalSystem::SparsityPattern::is_set_global_diagonal()
{
  return force_global_diagonal_;
}

void LocalSystem::SparsityPattern::resize(unsigned int nrows, unsigned int ncols)
{
    pattern.zeros(nrows,ncols);
}

const arma::umat & LocalSystem::SparsityPattern::get_matrix()
{
    return pattern;
}


