
#include "local_system.hh"

#include <armadillo>
#include "system/sys_vector.hh"

LocalSystem::LocalSystem(unsigned int nrows, unsigned int ncols)
: matrix(nrows, ncols), rhs(nrows)
{
    row_dofs.resize(matrix.n_rows);
    col_dofs.resize(matrix.n_cols);
    reset();
}


void LocalSystem::reset()
{
    // zeros in local system
    matrix.zeros();
    rhs.zeros();
    // drop all dirichlet values
    global_solution_dofs.clear();
    solution.clear();
    preferred_diag_values.clear();
    solution_not_set = true;
    
    // reset global degrees of freedom vectors
    std::fill(row_dofs.begin(), row_dofs.end(), 0);
    std::fill(col_dofs.begin(), col_dofs.end(), 0);
}


void LocalSystem::set_solution(unsigned int global_row, double solution_val, double diag_val)
{
//     ASSERT_DBG(loc_row < matrix.n_rows);
    global_solution_dofs.push_back(global_row);
    solution.push_back(solution_val);
    preferred_diag_values.push_back(diag_val);
    solution_not_set = false;
}


void LocalSystem::eliminate_solution()
{
    // skip diagonal fix
    if(solution_not_set) return;
    
//     DBGCOUT("fix_diagonal\n");
    
    arma::mat tmp_mat = matrix;
    arma::vec tmp_rhs = rhs;
    bool eliminate_row = false,
         eliminate_col = false;
    
    unsigned int i, l_row, l_col;
    
    // eliminate rows
    for(auto& sol_dof : global_solution_dofs)
        for(l_row = 0; l_row < matrix.n_rows; l_row++)
            if (row_dofs[l_row] == sol_dof) {
//                 DBGVAR(l_row);
                eliminate_row = true;
                tmp_rhs(l_row) = 0.0;
                tmp_mat.row(l_row).zeros();
            }
    
    // eliminate columns
    for(i = 0; i < global_solution_dofs.size(); i++){
        for(l_col = 0; l_col < matrix.n_cols; l_col++)
            if (col_dofs[l_col] == global_solution_dofs[i]) {
//                 DBGVAR(l_col);
                eliminate_col = true;
                tmp_rhs -= solution[i] * tmp_mat.col(l_col);
                tmp_mat.col(l_col).zeros();
            }
    }
    
    // correction of dirichlet diagonal entry
    // if both true, then there is eliminated diagonal entry
    if(eliminate_row && eliminate_col){
        unsigned int j, sol_dof;
        for(i=0; i < global_solution_dofs.size(); i++){
            sol_dof = global_solution_dofs[i];
            for(j = 0; j < matrix.n_rows; j++)  // find local row index of the solution
                if(sol_dof == row_dofs[j]){
                    l_row = j;
                    break;
                }
//             DBGVAR(sol_dof);
//             DBGVAR(l_row);
            for(l_col = 0; l_col < matrix.n_cols; l_col++){
//                 DBGVAR(col_dofs[l_col]);
                if (row_dofs[l_row] == col_dofs[l_col]){ // look for global diagonal entry
            
//                     DBGVAR(l_col);
                    // if preferred value is not set, then try using matrix value
                    double new_diagonal = matrix(l_row, l_col);
            
                    if(preferred_diag_values[i] !=0)    // if preferred value is set
                        new_diagonal = preferred_diag_values[i];
                    else if(new_diagonal == 0)      // if an assembled value is not available
                        new_diagonal = 1.0;
                    
//                     DBGVAR(new_diagonal);
                    
//                     double new_diagonal = fabs(matrix(sol_row,col));
//                     if (new_diagonal == 0.0) {
//                         if (matrix.is_square()) {
//                             new_diagonal = arma::sum( abs(matrix.diag())) / matrix.n_rows;
//                         } else {
//                             new_diagonal = arma::accu( abs(matrix) ) / matrix.n_elem;
//                         }
//                     }
                    tmp_mat(l_row,l_col) = new_diagonal;
                    tmp_rhs(l_row) = new_diagonal * solution[i];
                }
            }
        }
    }
    
    matrix = tmp_mat;
    rhs = tmp_rhs;
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