/*
 * linsys_test.cpp
 *
 *  Created on: Jan 16, 2014
 *      Author: jb
 */




#define TEST_USE_PETSC

#include "flow_gtest_mpi.hh"
#include "la/linsys.hh"
#include <armadillo>
#include "mpi.h"


class SetValues : public testing::Test, public LinSys {
public:
        /// LinSys implementation
	SetValues() : LinSys( new Distribution(10, MPI_COMM_WORLD) ) {
            matrix_.zeros();
            rhs_.zeros();
	}

	double get_solution_precision() override {
		return 0.0;
	}

	~SetValues() {
		delete rows_ds_;
	}

	void set_tolerances(double  r_tol, double a_tol, unsigned int max_it) override {
	}


	void mat_set_values(int nrow, int *rows, int ncol, int *cols,double *vals) {
		//cout << endl << "Matrix:" << endl;
		for(int i =0; i<nrow; i++) {
			for(int j=0; j<ncol; j++) {
                                this->matrix_(rows[i], cols[j])+=vals[i*ncol+j];
				//cout << "( " << rows[i] << ", " << cols[j] << ", " << vals[i*ncol+j] << ")" << endl;
			}
		}
	}

	void rhs_set_values(int nrow,int *rows,double *vals) {
		//cout << endl << "RHS:" << endl;
		for(int i =0; i<nrow; i++) {
                        this->rhs_( rows[i] )+=vals[i];
			//cout << "( " << rows[i] << ", " << vals[i] << ")" <<endl;
		}

	}

	void finish_assembly() {}

	void apply_constrains(double) {}

	int solve() {
		return 0;
	}
	
	
	/// Test methods
	/**
         * Set system size and generate random dirichlet conditions.
         */
	void set_size(unsigned int size) {
            full_matrix_=matrix_=arma::zeros(size, size);
            full_rhs_=rhs_=arma::zeros(size);
            
            //dirichlet_=arma::ones(size)-2*(arma::randu<arma::vec>(size)<0.2);
            dirichlet_rows_ = arma::find(arma::randu<arma::vec>(size)<0.2);
            dirichlet_.ones(size);
            dirichlet_.elem(dirichlet_rows_ )*=-1;
            dirichlet_[0]=1;
            
            dirichlet_values_=arma::randu<arma::vec>(size);
            dirichlet_rows_=arma::find(dirichlet_ < 0);
            non_dirichlet_rows_=arma::find(dirichlet_ > 0);
	}
	
	/**
         * Add a random local matrix and rhs spanning over given rows and columns.
         * 
         */
	void add(arma::uvec rows, arma::uvec cols) {
          
            arma::mat loc_mat=arma::randu<arma::mat>(rows.size(), cols.size());
            arma::vec loc_rhs=arma::randu<arma::vec>(rows.size());
            // apply to full system
            full_matrix_.submat(rows, cols)+=loc_mat;
            full_rhs_.elem(rows)+=loc_rhs;
            // apply to fixture system
            arma::vec row_sol=dirichlet_values_.elem(rows);
            arma::vec col_sol=dirichlet_values_.elem(cols);
            
            
            auto i_rows=arma::conv_to<std::vector<int> >::from(
                          arma::conv_to<arma::ivec>::from(rows)%dirichlet_.elem(rows));
            auto i_cols=arma::conv_to<std::vector<int> >::from(
                          arma::conv_to<arma::ivec>::from(cols)%dirichlet_.elem(cols));
            
            //cout << "i_rows\n" << arma::ivec(i_rows);
            //cout << "i_cols\n" << arma::ivec(i_cols);
            this->set_values(i_rows, i_cols, loc_mat, loc_rhs, row_sol, col_sol);            
            
            
            // check consistency
            double eps=4*arma::datum::eps;
            // zero dirichlet rows and cols
            //cout << "Dirich rows:\n" << dirichlet_rows_;
            //cout << "matrix_:\n" << matrix_;
            EXPECT_TRUE( arma::norm(matrix_.submat(dirichlet_rows_, non_dirichlet_rows_), "inf") < eps);
            
            EXPECT_TRUE( arma::norm(matrix_.submat(non_dirichlet_rows_, dirichlet_rows_), "inf") < eps);
            
            auto dirich_submat = matrix_.submat(dirichlet_rows_, dirichlet_rows_);
            EXPECT_TRUE( arma::norm( dirich_submat - arma::diagmat(dirich_submat), "inf") < eps );
            
            /*
            // full check
            arma::mat reduced_matrix_=full_matrix_;
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
	
	double compute_residual() {return 0;}

	
	arma::uvec non_dirichlet_rows_;
	arma::uvec dirichlet_rows_;
	arma::ivec dirichlet_;
        arma::vec dirichlet_values_;
        
	arma::mat matrix_;
        arma::vec rhs_;
        
        arma::mat full_matrix_;
        arma::vec full_rhs_;
};





TEST_F(SetValues, dirichlet) {
        
    for(unsigned int i=0; i<100; i++) {
        this->set_size(6);
        this->add( {0,1,2}, {0,1,2} );
        this->add( {0,1}, {3,4,5} );
        this->add( {4,5}, {0,1} );
        this->add( {0,2,3}, {0,2,3} );
        this->add( {3,4}, {3,4,5} );
        this->add( {5}, {0,2,4,5} );
        this->add( {1,3,4}, {4,5} );
        this->add( {0,3}, {4,5,} );
    }     
};
