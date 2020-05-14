/*
 * local_system_test.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: pe
 */



#include "flow_gtest.hh"
#include "la/local_system.hh"
#include <armadillo>
#include "arma_expect.hh"
#include "system/system.hh"
using namespace std;


/**
 * Simple example with rectangular local system (LS),
 * where the global diagonal is not aligned with the LS diagonal.
 * Test uses add_value(row,col,mat,rhs) function.
 * LocalSystem is used twice - once with preferred diagonal entries.
 * One dirichletBC condition has diagonal entry inside.
 * The second condition has diagonal entry outside the LS (makes zero row at the end).
 * 
 * local system in global:
 *    7  8  9  10  11
 * 10
 * 11                    dirBC s1
 * 12
 * 13                    dirBC s2
 */

void fill_ls(LocalSystem& ls){
    
    const unsigned int m = ls.get_matrix().n_rows,
                       n = ls.get_matrix().n_cols,
                       moff = 10,
                       noff = 7;
    // set global dofs
    for(unsigned int i=0; i < m; i++) ls.row_dofs[i] = moff + i;
    for(unsigned int i=0; i < n; i++) ls.col_dofs[i] = noff + i;
    
    //fill rows
    for(unsigned int i=0; i < n; i++){
        ls.add_value(0, i, 10*i, 50);
        ls.add_value(1, i, 10*(i+2), 60);
        ls.add_value(2, i, 10*(i+4), 70);
        ls.add_value(3, i, 10*(i+5));
        ls.add_value(3, 80);
    }
}

TEST(la, simple_local_system_solution) {

    const unsigned int m = 4,
                       n = 5,
                       moff = 10,
                       noff = 7;
    const double 
        s1 = 5,
        s2 = 6,
        d1 = 100,
        d2 = 200;
    arma::mat res_mat(m,n);
    res_mat.row(0) = arma::rowvec({0.0, 10, 20, 30, 0});
    res_mat.row(1) = arma::rowvec({0.0, 0, 0, 0, d1});
    res_mat.row(2) = arma::rowvec({40.0, 50, 60, 70, 0});
    res_mat.row(3) = arma::rowvec({0.0, 0, 0, 0, 0});
    arma::vec res_rhs = {n*50.0 - 40*s1,
                         d1*5.0,
                         n*70.0 - 80*s1,
                         0};
    
    LocalSystem ls(m, n);
    
    // set sparsity pattern
    arma::umat sp;
    sp.ones(m,n);
    ls.set_sparsity(sp);
    
    // set solution with preferred diagonal entries
    //ls.set_solution(11, s1, d1);
    ls.set_solution_row(1, s1, d1);
    ls.set_solution_col(4, s1);
    //ls.set_solution(13, s2, d2);
    ls.set_solution_row(3, s2, d2);
    //ls.set_solution_col(6, s1); // out of matrix
    
    fill_ls(ls);
//     cout << "matrix:\n" << ls.get_matrix();
//     cout << "rhs:\n" << ls.get_rhs();
    
    ls.eliminate_solution();
    
//     cout << "matrix:\n" << ls.get_matrix();
//     cout << "res_mat:\n" << res_mat;
//     cout << "rhs:\n" << ls.get_rhs();
//     cout << "res_rhs:\n" << res_rhs;
    
    EXPECT_ARMA_EQ(res_mat, ls.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, ls.get_rhs());
    
    //change and do not set the preferred diagonal value
    ls.reset();
    ls.set_solution_row(1, s1);
    ls.set_solution_col(4, s1);
    ls.set_solution_row(3, s2);

    fill_ls(ls);
    ls.eliminate_solution();
    
    res_mat(1,4) = 60;
    res_rhs(1) = 60*s1;
//     cout << "matrix:\n" << ls.get_matrix();
//     cout << "res_mat:\n" << res_mat;
    
    EXPECT_ARMA_EQ(res_mat, ls.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, ls.get_rhs());
    
    
    //same as before, but add solution, that corresponds only with column 8 (should be eliminated)
    ls.reset();
    ls.set_solution_row(1, s1);
    ls.set_solution_col(4, s1);

    ls.set_solution_row(3, s2);

    ls.set_solution_col(1, s2);

    fill_ls(ls);
    ls.eliminate_solution();
    
    res_mat.col(1).zeros();
    res_rhs(0) -= 10*s2;
    res_rhs(2) -= 50*s2;
//     cout << "matrix:\n" << ls.get_matrix();
//     cout << "res_mat:\n" << res_mat;
    
    EXPECT_ARMA_EQ(res_mat, ls.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, ls.get_rhs());
}




TEST(la, simple_local_system_no_solution) {

    const unsigned int m = 2,
                       n = 3;
    arma::mat res_mat(m,n);
    res_mat.row(0) = arma::rowvec({1.0, 2, 3});
    res_mat.row(1) = arma::rowvec({4.0, 5, 6});
    arma::vec res_rhs = {10, 11};
    
    LocalSystem ls(m, n);
    
    // set sparsity pattern
    arma::umat sp;
    sp.ones(m,n);
    ls.set_sparsity(sp);
    
    ls.add_value(0, 0, 1, 10);
    ls.add_value(0, 1, 2, 0);
    ls.add_value(0, 2, 3, 0);
    ls.add_value(1, 0, 4, 11);
    ls.add_value(1, 1, 5, 0);
    ls.add_value(1, 2, 6, 0);
    
    ls.eliminate_solution();
    
    EXPECT_ARMA_EQ(res_mat, ls.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, ls.get_rhs());
}


TEST(la, sparse_local_system) {

    const unsigned int m = 4,
                       n = 4;
    // set sparsity pattern
    arma::umat sp(m,n);
    sp.zeros();
    sp.diag() = arma::uvec({1,2,3,4});
    sp.diag(2) = arma::uvec({5,6});
    sp.diag(-3) = arma::uvec({7});
    
    //resulting system
    arma::mat res_mat(m,n);
    for(unsigned int i=0; i<m;i++)
        for(unsigned int j=0; j<n;j++)
            res_mat(i,j) = sp(i,j);
    arma::vec res_rhs = {10, 11, 12, 13};
    
    LocalSystem ls(m, n);
    
    // set sparsity pattern
    ls.set_sparsity(sp);
    
    EXPECT_ASSERT_DEATH(ls.add_value(0, 1, 1.0), "Violation of sparsity pattern.");
    EXPECT_ASSERT_DEATH(ls.add_value(0, 3, 1.0), "Violation of sparsity pattern.");
    EXPECT_ASSERT_DEATH(ls.add_value(2, 0, 1.0, 1.0), "Violation of sparsity pattern.");
    EXPECT_ASSERT_DEATH(ls.add_value(3, 1, 1.0, 1.0), "Violation of sparsity pattern.");
    ls.add_value(0, 0, 1, 10);
    ls.add_value(1, 1, 2, 11);
    ls.add_value(2, 2, 3, 12);
    ls.add_value(3, 3, 4, 13);
    ls.add_value(0, 2, 5);
    ls.add_value(1, 3, 6);
    ls.add_value(3, 0, 7);
    
    ls.eliminate_solution();
    
    EXPECT_ARMA_EQ(res_mat, ls.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, ls.get_rhs());
}











class SetValues : public testing::Test, public LocalSystem {
public:
    
    /// fixed size of the local system
    static const unsigned int size = 9;
    
    SetValues() : LocalSystem(size,size) {
    }

    ~SetValues() {
    }

    /// Test methods
    /**
     * Set system size and generate random dirichlet conditions.
     */
    void restart() {
            full_matrix_ = matrix = arma::zeros(size, size);
            full_rhs_    = rhs    = arma::zeros(size);
            
            // set sparsity pattern
            arma::umat sp;
            sp.ones(size, size);
            set_sparsity(sp);
            
            dirichlet_rows_ = arma::find(arma::randu<arma::vec>(size)<0.2);
            dirichlet_.ones(size);
            dirichlet_.elem(dirichlet_rows_ ) *= -1;
            dirichlet_[0]=1;
            
            dirichlet_values_=arma::randu<arma::vec>(size);
            dirichlet_rows_=arma::find(dirichlet_ < 0);
            non_dirichlet_rows_=arma::find(dirichlet_ > 0);
            
            preferred_diag_ = arma::randu<arma::vec>(dirichlet_rows_.size());
            
//             cout << "dir_rows\n" << dirichlet_rows_;
//             cout << "dir\n" << dirichlet_;
            
            this->reset();
            for(unsigned int i=0; i < size; i++){
                row_dofs[i] = i;
                col_dofs[i] = i;
            }
            
            // set dirichlet BC
            for(unsigned int i=0; i < dirichlet_rows_.n_elem; i++){
//                 cout << "set solution: " << dirichlet_rows_(i) << "  " << dirichlet_values_(dirichlet_rows_(i)) << endl;
                if(preferred_flag_)
                    this->set_solution(dirichlet_rows_(i), dirichlet_values_(dirichlet_rows_(i)), preferred_diag_(i));
                else
                    this->set_solution(dirichlet_rows_(i), dirichlet_values_(dirichlet_rows_(i)));
            }   
    }
    
    /** 
     * Generate random local matrix and rhs, which is then added into the full_matrix and full_rhs
     * according to the given rows and columns indices.
     */
    void edit_full_matrix(arma::uvec rows, arma::uvec cols, arma::mat& loc_mat, arma::vec& loc_rhs, bool add){
        loc_mat=arma::randu<arma::mat>(rows.size(), cols.size());
        loc_rhs=arma::randu<arma::vec>(rows.size());
        // apply to full system
        if(add){
            full_matrix_.submat(rows, cols)+=loc_mat;
            full_rhs_.elem(rows)+=loc_rhs;
        }
        else {
            full_matrix_.submat(rows, cols)=loc_mat;
            full_rhs_.elem(rows)=loc_rhs;
        }
            
//         cout << "full_matrix\n" << full_matrix_;
//         cout << "full_rhs\n" << full_rhs_;
    }
    
    /** 
     * Performs checking of the dirichlet elimination.
     */
    void check_result(bool full = false)
    {
        
        // print results:
//         cout << "Dirich rows:\n" << dirichlet_rows_;
//         cout << "Dirich values:\n" << dirichlet_values_;
//         cout << "preferred_diag_:\n" << preferred_diag_;
//         cout << "matrix_:\n" << matrix;
//         cout << "full_matrix_:\n" << full_matrix_;
//         cout << "rhs:\n" << rhs;
//         cout << "full_rhs:\n" << full_rhs_;
//         cout << "diag_value:\n";
//         for(unsigned int i=0; i < diag_values.size(); i++) cout << diag_values[i] << "\n";
        
        // check consistency
        double eps = 8*arma::datum::eps;
        
        // extracts the rows and columns which are supposed to be eliminated by dirichlet BC
        EXPECT_LT( arma::norm(matrix.submat(dirichlet_rows_, non_dirichlet_rows_), "inf"), eps);    
        EXPECT_LT( arma::norm(matrix.submat(non_dirichlet_rows_, dirichlet_rows_), "inf"), eps);
        
        // extracts matrix of dirichletBC rows and checks that it is diagonal
        arma::mat dirich_submat = matrix.submat(dirichlet_rows_, dirichlet_rows_);
        EXPECT_LT( arma::norm( dirich_submat - arma::diagmat(dirich_submat), "inf"), eps );
        
        // check Dirichlet BC diagonal entries and rhs:
        if(preferred_flag_){  //if preferred values are set
            arma::vec d_rhs = rhs.elem(dirichlet_rows_);
            for(unsigned int i=0; i < d_rhs.size(); i++){
                EXPECT_EQ(dirich_submat(i,i), preferred_diag_(i));
                EXPECT_EQ(d_rhs(i), preferred_diag_(i) * dirichlet_values_(dirichlet_rows_(i)));
            }
        }
        else {  //if not, then diag values are summed or is equal 1.0
            arma::vec d_rhs = rhs.elem(dirichlet_rows_);
            for(unsigned int i=0; i < d_rhs.size(); i++){
                EXPECT_EQ(d_rhs(i), dirich_submat(i,i) * dirichlet_values_(dirichlet_rows_(i)));
//                 if(diag_values[dirichlet_rows_(i)] == 0){
// //                     cout << "DIAG VALUES 1.0\n";
//                     EXPECT_EQ(dirich_submat(i,i), 1.0);
//                     EXPECT_EQ(d_rhs(i), dirichlet_values_[dirichlet_rows_(i)]);
//                 }
//                 else {
//                     cout << "DIAG VALUES\n";
//                     EXPECT_EQ(dirich_submat(i,i), diag_values[dirichlet_rows_(i)]);
//                     EXPECT_EQ(d_rhs(i), diag_values[dirichlet_rows_(i)] * dirichlet_values_(dirichlet_rows_(i)));
//                 }
            }
        }
        
        
        if(full) { // full check
            // dirich_submat is checked at this position
            
            arma::uvec vec_span(size); for(unsigned int i=0; i < size; i++) vec_span(i) = i;
            arma::mat reduced_mat = full_matrix_;
            arma::mat dirich_cols = reduced_mat.submat(vec_span, dirichlet_rows_);
            
            // eliminate dirichletBC from rhs:
            arma::vec reduced_rhs = full_rhs_;
            for(unsigned int i=0; i < dirich_cols.n_cols; i++) 
                reduced_rhs = reduced_rhs - dirichlet_values_(dirichlet_rows_(i)) * dirich_cols.col(i);
            
//             cout << "dirich_cols:\n" << dirich_cols;
//             cout << "rhs:\n" << rhs;
//             cout << "reduced rhs:\n" << reduced_rhs;
            EXPECT_LT( arma::norm(rhs.elem(non_dirichlet_rows_) - reduced_rhs.elem(non_dirichlet_rows_), "inf"), eps );
            
            reduced_rhs.elem(dirichlet_rows_) = dirich_submat * dirichlet_values_.elem(dirichlet_rows_);
            EXPECT_LT( arma::norm(rhs - reduced_rhs, "inf"), eps );
            
            reduced_mat.submat(dirichlet_rows_, vec_span).zeros();
            reduced_mat.submat(vec_span, dirichlet_rows_).zeros();
            reduced_mat.submat(dirichlet_rows_, dirichlet_rows_) = dirich_submat;
            EXPECT_LT( arma::norm(matrix - reduced_mat, "inf"), eps );
        }
    }
    
    
    /**
     * Adds a random local matrix and rhs spanning over given rows and columns.
     * Add value one by one.
     */
    void add_value_single(arma::uvec rows, arma::uvec cols) {
        
        arma::mat loc_mat;
        arma::vec loc_rhs;
        edit_full_matrix(rows,cols, loc_mat, loc_rhs, true);
        
        // set entries
        for(unsigned int i=0; i < rows.n_elem; i++){
            this->add_value(rows(i), 0, 0.0, loc_rhs(i));
            for(unsigned int j=0; j < cols.n_elem; j++){
                this->add_value(rows(i), cols(j), loc_mat(i,j), 0.0);
            }
        }
    }

    void set_preferred()
    { preferred_flag_ = true; }
    
    void unset_preferred()
    { preferred_flag_ = false; }
    
    arma::uvec non_dirichlet_rows_; ///< indices of rows with no dirichletBC
    arma::uvec dirichlet_rows_;     ///< indices of rows with dirichletBC
    arma::ivec dirichlet_;          ///< indices of rows, dirichletBC marked by negative index (for original function by JB)
    arma::vec dirichlet_values_;    ///< values of DoFs set by dirichletBC
    arma::vec preferred_diag_;      ///< values that are to be set at the diagonal entry of rows with dirichletBC
    bool preferred_flag_;
        
    arma::mat full_matrix_;         ///< full matrix with no dirichletBC elimination
    arma::vec full_rhs_;            ///< full rhs with no dirichletBC elimination
};


TEST_F(SetValues, add_value_single) {
    
    unset_preferred();
    unsigned int n = 100;
    for(unsigned int i=0; i<n; i++) 
    {
//         cout << "############################################################   " << i << endl;
        restart();
        add_value_single( {0,1,2}, {0,1,2} );
        add_value_single( {0,1}, {3,4,5} );
        add_value_single( {4,5}, {0,1} );
        add_value_single( {0,2,3}, {0,2,3} );
        add_value_single( {3,4}, {3,4,5} );
        add_value_single( {5}, {0,2,4,5} );
        add_value_single( {1,3,4}, {4,5} );
        add_value_single( {0,3}, {4,5} );
        
        this->eliminate_solution();
        check_result(true);
    }
    
    set_preferred();
    for(unsigned int i=0; i<n; i++) 
    {
//         cout << "############################################################   " << i << endl;
        restart();
        add_value_single( {0,1,2}, {0,1,2} );
        add_value_single( {0,1}, {3,4,5} );
        add_value_single( {4,5}, {0,1} );
        add_value_single( {0,2,3}, {0,2,3} );
        add_value_single( {3,4}, {3,4,5} );
        add_value_single( {5}, {0,2,4,5} );
        add_value_single( {1,3,4}, {4,5} );
        add_value_single( {0,3}, {4,5} );
        
        this->eliminate_solution();
        check_result(true);
    }  
};


TEST(la, schur_complement) {

    const unsigned int m = 5,
                       n = 5;
    
    arma::mat M = {{1, 1, -1, 1, 2}, {1, 2, 1, 2, 0}, {2, -1, 1, 3, 1},
                   {1, 2, 3, 4, 1}, {2, 0, 1, 1, 2}};
//     M.print();
    arma::vec rhs = {1, -2, 1, 2, -1};
    
    // set sparsity pattern
    arma::umat sp(m,n);
    sp.ones();
    
    LocalSystem ls(m, n);
    
    // set sparsity pattern
    ls.set_sparsity(sp);
    
    ls.set_matrix(M);
    ls.set_rhs(rhs);
    
//     ls.eliminate_solution();
    LocalSystem schur;
    ls.compute_schur_complement(3, schur);
    
//     schur.get_matrix().print();
//     schur.get_rhs().print();
    
    arma::mat res_mat = {{1+1./9, 3}, {-2-1./9, 1}};
    arma::vec res_rhs = {6+2./9, -1-2./9};
    EXPECT_ARMA_EQ(res_mat, schur.get_matrix());
    EXPECT_ARMA_EQ(res_rhs, schur.get_rhs());
    
    
    arma::vec schur_sol = arma::solve(res_mat, res_rhs);
    arma::vec res_sol = arma::solve(M, rhs);
    
    arma::vec reconstructed_solution;
    ls.reconstruct_solution_schur(3, schur_sol, reconstructed_solution);
//     reconstructed_solution.print();
    EXPECT_ARMA_EQ(res_sol.subvec(0,2), reconstructed_solution);
}
