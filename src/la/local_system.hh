
#ifndef LOCAL_SYSTEM_HH_
#define LOCAL_SYSTEM_HH_

#include <armadillo>
#include "system/global_defs.h"

class LocalSystem
    {
    public:
        
        LocalSystem(unsigned int nrows, unsigned int ncols)
        : matrix(arma::zeros(nrows, ncols)), rhs(arma::zeros(nrows,1))
        {}
        
        void reset(){
            // zeros in local system
            matrix.zeros();
            rhs.zeros();
            // drop all dirichlet values
            loc_dirichlet_edge.clear();
            loc_dirichlet_val.clear();
            // reset global degrees of freedom vectors
            row_dofs.resize(row_dofs.size(),0);
            col_dofs.resize(col_dofs.size(),0);
        }
        
        arma::mat matrix;
        arma::vec rhs;
        
        std::vector<unsigned int> loc_dirichlet_edge;   //rename to dirichlet_indicator
        std::vector<double> loc_dirichlet_val;
        
        std::vector<int> row_dofs;
        std::vector<int> col_dofs;
        
        
        
        
        /**
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
        void set_values(std::vector<int> &row_dofs, std::vector<int> &col_dofs,
                        const arma::mat &loc_matrix, const arma::vec &loc_rhs,
                        const arma::vec &row_solution, const arma::vec &col_solution)

        {
            arma::mat tmp = loc_matrix;
            arma::vec tmp_rhs = loc_rhs;
            bool negative_row = false;
            bool negative_col = false;

//             DBGCOUT(<< "tmp\n" << tmp);
//             DBGCOUT("lrow\n");
            for(unsigned int l_row = 0; l_row < row_dofs.size(); l_row++)
                if (row_dofs[l_row] < 0) {
                    tmp_rhs(l_row)=0.0;
                    tmp.row(l_row).zeros();
                    negative_row=true;
                }

//             DBGCOUT("lcol\n");
            for(unsigned int l_col = 0; l_col < col_dofs.size(); l_col++)
                if (col_dofs[l_col] < 0) {
                    tmp_rhs -= loc_matrix.col(l_col) * col_solution[l_col];
                    tmp.col(l_col).zeros();
                    negative_col=true;
                }
                
//             DBGCOUT("main\n");
                
            if (negative_row && negative_col) {
                // look for diagonal entry
                for(unsigned int l_row = 0; l_row < row_dofs.size(); l_row++)
                    if (row_dofs[l_row] < 0)
                        for(unsigned int l_col = 0; l_col < col_dofs.size(); l_col++)
                            if (col_dofs[l_col] < 0 && row_dofs[l_row] == col_dofs[l_col]) {
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
                for(int &row : row_dofs) row=abs(row);

            if (negative_col)
                for(int &col : col_dofs) col=abs(col);

            
//             DBGCOUT( << "row_dofs:\n" << arma::conv_to<arma::uvec>::from(row_dofs));
//             DBGCOUT( << "matrix:\n" << matrix);
//             DBGCOUT("set mat and rhs\n");
//             DBGCOUT(<< "tmp\n" << tmp);
            matrix.submat(arma::conv_to<arma::uvec>::from(row_dofs), arma::conv_to<arma::uvec>::from(col_dofs)) = tmp;
            rhs.elem(arma::conv_to<arma::uvec>::from(row_dofs)) = tmp_rhs;
        }
    };
    
#endif // LOCAL_SYSTEM_HH_