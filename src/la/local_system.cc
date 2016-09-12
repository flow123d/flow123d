
#include "local_system.hh"

#include <armadillo>
#include "system/global_defs.h"

LocalSystem::LocalSystem(unsigned int nrows, unsigned int ncols)
: matrix(nrows, ncols), rhs(nrows)
{
    reset();
}


void LocalSystem::reset()
{
    // zeros in local system
    matrix.zeros();
    rhs.zeros();
    // drop all dirichlet values
    loc_solution_rows.clear();
    loc_solution.clear();
    // reset global degrees of freedom vectors
    row_dofs.resize(matrix.n_rows,0);
    col_dofs.resize(matrix.n_cols,0);
}


void LocalSystem::set_solution(unsigned int loc_row, double solution)
{
    ASSERT_DBG(loc_row < matrix.n_rows);
    loc_solution_rows.push_back(loc_row);
    loc_solution.push_back(solution);
}


void LocalSystem::fix_diagonal()
{
    // correction of dirichlet diagonal entry
    for(auto& sol_row : loc_solution_rows)
        for(unsigned int col = 0; col < matrix.n_cols; col++) {
//                     DBGCOUT(<< "row " << row_dofs[sol_row] << "  col " << col_dofs[col] << "\n");
            // look for diagonal entry
            if (row_dofs[sol_row] == col_dofs[col]) {
//                                 DBGCOUT(<< "lrow=" << l_row << " lcol=" << l_col << "\n");
                double new_diagonal = fabs(matrix(sol_row,col));
                if (new_diagonal == 0.0) {
                    if (matrix.is_square()) {
                        new_diagonal = arma::sum( abs(matrix.diag())) / matrix.n_rows;
                    } else {
                        new_diagonal = arma::accu( abs(matrix) ) / matrix.n_elem;
                    }
                }
                matrix(sol_row,col) = new_diagonal;
                rhs(sol_row) = new_diagonal * loc_solution[sol_row];
            }
        }
}
        
void LocalSystem::set_value(unsigned int row, unsigned int col, double mat_val, double rhs_val)
{
//             DBGCOUT(<< "row " << row_dof << "  col " << col_dof << "\n");
    ASSERT_DBG(row < matrix.n_rows);
    ASSERT_DBG(col < matrix.n_cols);
    
    bool eliminate_row = false;
    
    double tmp_mat = mat_val;
    double tmp_rhs = rhs_val;
    
    for(auto& sol : loc_solution_rows) {
        if(sol == col){
            tmp_mat = 0.0;
            tmp_rhs -= mat_val * loc_solution[col];
        }
        if(sol == row){ eliminate_row = true; }
    }
    
    if(! eliminate_row){
        matrix(row, col) = tmp_mat;
        rhs(row) = tmp_rhs;
    }
}


void LocalSystem::set_values(std::vector< unsigned int >& rows, std::vector< unsigned int >& cols,
                             const arma::mat& loc_matrix, const arma::vec& loc_rhs)
{
    ASSERT_DBG(loc_matrix.n_rows <= matrix.n_rows);
    ASSERT_DBG(loc_matrix.n_cols <= matrix.n_cols);
    
    bool eliminate_row = false;
    
    arma::mat tmp_mat = loc_matrix;
    arma::vec tmp_rhs = loc_rhs;
    
//     DBGCOUT("lrow\n");
    for(auto& sol : loc_solution_rows)
        for(unsigned int l_row = 0; l_row < rows.size(); l_row++)
            if (rows[l_row] == sol) {
                tmp_rhs(l_row) = 0.0;
                tmp_mat.row(l_row).zeros();
            }
    
//     DBGCOUT("lcol\n");
    for(auto& sol : loc_solution_rows)
        for(unsigned int l_col = 0; l_col < cols.size(); l_col++)
            if (cols[l_col] == sol) {
                tmp_rhs -= loc_matrix.col(l_col) * loc_solution[sol];
                tmp_mat.col(l_col).zeros();
            }
    
//     DBGCOUT("set mat and rhs\n");
//     DBGCOUT(<< "tmp\n" << tmp_mat);
    matrix.submat(arma::conv_to<arma::uvec>::from(rows), arma::conv_to<arma::uvec>::from(cols)) = tmp_mat;
    rhs.elem(arma::conv_to<arma::uvec>::from(rows)) = tmp_rhs;
}


void LocalSystem::set_values(std::vector< int >& rows, std::vector< int >& cols,
                             const arma::mat& loc_matrix, const arma::vec& loc_rhs,
                             const arma::vec& row_solution, const arma::vec& col_solution)
{
    ASSERT_DBG(loc_matrix.n_rows <= matrix.n_rows);
    ASSERT_DBG(loc_matrix.n_cols <= matrix.n_cols);
    
    arma::mat tmp = loc_matrix;
    arma::vec tmp_rhs = loc_rhs;
    bool negative_row = false;
    bool negative_col = false;

//             DBGCOUT(<< "tmp\n" << tmp);
//             DBGCOUT("lrow\n");
    for(unsigned int l_row = 0; l_row < rows.size(); l_row++)
        if (rows[l_row] < 0) {
            tmp_rhs(l_row)=0.0;
            tmp.row(l_row).zeros();
            negative_row=true;
        }

//             DBGCOUT("lcol\n");
    for(unsigned int l_col = 0; l_col < cols.size(); l_col++)
        if (cols[l_col] < 0) {
            tmp_rhs -= loc_matrix.col(l_col) * col_solution[l_col];
            tmp.col(l_col).zeros();
            negative_col=true;
        }
        
//             DBGCOUT("main\n");
        
    if (negative_row && negative_col) {
        // look for diagonal entry
        for(unsigned int l_row = 0; l_row < rows.size(); l_row++)
            if (rows[l_row] < 0)
                for(unsigned int l_col = 0; l_col < cols.size(); l_col++)
                    if (cols[l_col] < 0 && row_dofs[-rows[l_row]] == col_dofs[-cols[l_col]]) {
//                                 DBGCOUT(<< "lrow=" << l_row << " lcol=" << l_col << "\n");
                        double new_diagonal = fabs(loc_matrix.at(l_row,l_col));
                        if (new_diagonal == 0.0) {
                            if (loc_matrix.is_square()) {
                                new_diagonal = arma::sum( abs(loc_matrix.diag())) / loc_matrix.n_rows;
                            } else {
                                new_diagonal = arma::accu( abs(loc_matrix) ) / loc_matrix.n_elem;
                            }
                        }
//                                 tmp.at(l_col, l_row) = new_diagonal;
                        tmp.at(l_row, l_col) = new_diagonal;
                        tmp_rhs(l_row) = new_diagonal * row_solution[l_row];
                    }

    }

    if (negative_row)
        for(int &row : rows) row=abs(row);

    if (negative_col)
        for(int &col : cols) col=abs(col);

    
//             DBGCOUT( << "row_dofs:\n" << arma::conv_to<arma::uvec>::from(row_dofs));
//             DBGCOUT( << "matrix:\n" << matrix);
//             DBGCOUT("set mat and rhs\n");
//             DBGCOUT(<< "tmp\n" << tmp);
    matrix.submat(arma::conv_to<arma::uvec>::from(rows), arma::conv_to<arma::uvec>::from(cols)) = tmp;
    rhs.elem(arma::conv_to<arma::uvec>::from(rows)) = tmp_rhs;
}