/*
 * local_system_test.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: pe
 */



#include "flow_gtest.hh"
#include "la/local_system.hh"
#include <armadillo>

using namespace std;



class SetValues : public testing::Test, public LocalSystem {
public:
    
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
            
            dirichlet_rows_ = arma::find(arma::randu<arma::vec>(size)<0.2);
            dirichlet_.ones(size);
            dirichlet_.elem(dirichlet_rows_ ) *= -1;
            dirichlet_[0]=1;
            
            dirichlet_values_=arma::randu<arma::vec>(size);
            dirichlet_rows_=arma::find(dirichlet_ < 0);
            non_dirichlet_rows_=arma::find(dirichlet_ > 0);
            
//             cout << "dir_rows\n" << dirichlet_rows_;
//             cout << "dir\n" << dirichlet_;
            
            this->reset();
            for(unsigned int i=0; i < size; i++){
                row_dofs[i] = i;
                col_dofs[i] = i;
            }
            
    }
    
    void edit_full_matrix(arma::uvec rows, arma::uvec cols, arma::mat& loc_mat, arma::vec& loc_rhs){
        loc_mat=arma::randu<arma::mat>(rows.size(), cols.size());
        loc_rhs=arma::randu<arma::vec>(rows.size());
        // apply to full system
        full_matrix_.submat(rows, cols)+=loc_mat;
        full_rhs_.elem(rows)+=loc_rhs;
            
//         cout << "full_matrix\n" << full_matrix_;
//         cout << "full_rhs\n" << full_rhs_;
    }
    
    void check_result()
    {
        
        // zero dirichlet rows and cols
//         cout << "Dirich rows:\n" << dirichlet_rows_;
//         cout << "matrix_:\n" << matrix;
        
        // check consistency
        double eps=4*arma::datum::eps;
        
        EXPECT_TRUE( arma::norm(matrix.submat(dirichlet_rows_, non_dirichlet_rows_), "inf") < eps);
            
        EXPECT_TRUE( arma::norm(matrix.submat(non_dirichlet_rows_, dirichlet_rows_), "inf") < eps);
        
        auto dirich_submat = matrix.submat(dirichlet_rows_, dirichlet_rows_);
        EXPECT_TRUE( arma::norm( dirich_submat - arma::diagmat(dirich_submat), "inf") < eps );
        
        /*
        // full check
        arma::mat reduced_mat=full_matrix_;
        auto dirich_cols=reduced_matrix_.submat(arma::span::all, dirichlet_rows_);
        arma::vec reduced_rhs_=full_rhs_ - dirich_cols *dirichlet_values_;
        dirich_cols.zeros();
        
        reduced_matrix_.submat(dirichlet_rows_, arma::span::all).zeros();
        reduced_matrix_.submat(dirichlet_rows_, dirichlet_rows_) = dirich_submat;
        reduced_rhs_.subvec(dirichlet_rows_)=dirich_submat*dirichlet_values_;
        
        EXPECT_TRUE( arma::all( abs(matrix_ - reduced_matrix_)<eps ) );
        EXPECT_TRUE( arma::all( abs(rhs_ - reduced_rhs_)<eps ) );
        */
    }
    
    /**
         * Add a random local matrix and rhs spanning over given rows and columns.
         * 
         */
    void add(arma::uvec rows, arma::uvec cols) {
        
        arma::mat loc_mat;
        arma::vec loc_rhs;
        edit_full_matrix(rows,cols, loc_mat, loc_rhs);
        
        // apply to fixture system
        arma::vec row_sol=dirichlet_values_.elem(rows);
        arma::vec col_sol=dirichlet_values_.elem(cols);
        
        
        auto i_rows=arma::conv_to<std::vector<int> >::from(
                        arma::conv_to<arma::ivec>::from(rows)%dirichlet_.elem(rows));
        auto i_cols=arma::conv_to<std::vector<int> >::from(
                        arma::conv_to<arma::ivec>::from(cols)%dirichlet_.elem(cols));
        
//         cout << "i_rows\n" << arma::ivec(i_rows);
//         cout << "i_cols\n" << arma::ivec(i_cols);
        this->set_values(i_rows, i_cols, loc_mat, loc_rhs, row_sol, col_sol);            
        
        check_result();
    }
    
    /**
         * Add a random local matrix and rhs spanning over given rows and columns.
         * 
         */
    void add_single(arma::uvec rows, arma::uvec cols) {
        
        arma::mat loc_mat;
        arma::vec loc_rhs;
        edit_full_matrix(rows,cols, loc_mat, loc_rhs);
            
        // set dirichlet
        for(unsigned int i=0; i < dirichlet_rows_.n_elem; i++){
//             cout << "set solution: " << dirichlet_rows_(i) << "  " << dirichlet_values_(dirichlet_rows_(i)) << endl;
            this->set_solution(dirichlet_rows_(i), dirichlet_values_(dirichlet_rows_(i)));
        }
        
        // set entries
        for(unsigned int i=0; i < rows.n_elem; i++)
            for(unsigned int j=0; j < cols.n_elem; j++)
                this->set_value(rows(i), cols(j), loc_mat(i,j), loc_rhs(i));
        
        this->fix_diagonal();
        
        check_result();
    }
    
        /**
         * Add a random local matrix and rhs spanning over given rows and columns.
         * 
         */
    void add_subsystem(arma::uvec rows, arma::uvec cols) {
        
        arma::mat loc_mat;
        arma::vec loc_rhs;
        edit_full_matrix(rows,cols, loc_mat, loc_rhs);
        
//         auto i_rows=arma::conv_to<std::vector<int> >::from(
//                         arma::conv_to<arma::ivec>::from(rows)%dirichlet_.elem(rows));
//         auto i_cols=arma::conv_to<std::vector<int> >::from(
//                         arma::conv_to<arma::ivec>::from(cols)%dirichlet_.elem(cols));
        
        auto i_rows = arma::conv_to<std::vector<unsigned int>>::from(rows);
        auto i_cols = arma::conv_to<std::vector<unsigned int>>::from(cols);
//         cout << "i_rows\n" << arma::ivec(i_rows);
//         cout << "i_cols\n" << arma::ivec(i_cols);
        
        
        // set dirichlet
        for(unsigned int i=0; i < dirichlet_rows_.n_elem; i++){
//             cout << "set solution: " << dirichlet_rows_(i) << "  " << dirichlet_values_(dirichlet_rows_(i)) << endl;
            this->set_solution(dirichlet_rows_(i), dirichlet_values_(dirichlet_rows_(i)));
        }
        
        this->set_values(i_rows, i_cols, loc_mat, loc_rhs); 
        
        this->fix_diagonal();
        
        check_result();
    }

    
    arma::uvec non_dirichlet_rows_;
    arma::uvec dirichlet_rows_;
    arma::ivec dirichlet_;
    arma::vec dirichlet_values_;
        
    arma::mat full_matrix_;
    arma::vec full_rhs_;
};





TEST_F(SetValues, dirichlet_values) {
    
    unsigned int n = 100;
    for(unsigned int i=0; i<n; i++) {
//         cout << "############################################################   " << i << endl;
        restart();
        add( {0,1,2}, {0,1,2} );
        add( {0,1}, {3,4,5} );
        add( {4,5}, {0,1} );
        add( {0,2,3}, {0,2,3} );
        add( {3,4}, {3,4,5} );
        add( {5}, {0,2,4,5} );
        add( {1,3,4}, {4,5} );
        add( {0,3}, {4,5} );
    }     
};

TEST_F(SetValues, set_value) {
    
    unsigned int n = 100;
    for(unsigned int i=0; i<n; i++) {
//         cout << "############################################################   " << i << endl;
        restart();
        add_single( {0,1,2}, {0,1,2} );
        add_single( {0,1}, {3,4,5} );
        add_single( {4,5}, {0,1} );
        add_single( {0,2,3}, {0,2,3} );
        add_single( {3,4}, {3,4,5} );
        add_single( {5}, {0,2,4,5} );
        add_single( {1,3,4}, {4,5} );
        add_single( {0,3}, {4,5} );
    }     
}; 

TEST_F(SetValues, set_values) {
    
    unsigned int n = 100;
    for(unsigned int i=0; i<n; i++) {
//         cout << "############################################################   " << i << endl;
        restart();
        add_subsystem( {0,1,2}, {0,1,2} );
        add_subsystem( {0,1}, {3,4,5} );
        add_subsystem( {4,5}, {0,1} );
        add_subsystem( {0,2,3}, {0,2,3} );
        add_subsystem( {3,4}, {3,4,5} );
        add_subsystem( {5}, {0,2,4,5} );
        add_subsystem( {1,3,4}, {4,5} );
        add_subsystem( {0,3}, {4,5} );
    }     
};