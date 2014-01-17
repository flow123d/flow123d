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
	SetValues() : LinSys( new Distribution(10, MPI_COMM_WORLD) ) {
	}

	~SetValues() {
		delete rows_ds_;
	}

	void mat_set_values(int nrow,int *rows,int ncol,int *cols,double *vals) {
		cout << endl << "Matrix:" << endl;
		for(int i =0; i<nrow; i++) {
			for(int j=0; j<ncol; j++) {
				cout << "( " << rows[i] << ", " << cols[j] << ", " << vals[i*ncol+j] << ")" << endl;
			}
		}
	}

	void rhs_set_values(int nrow,int *rows,double *vals) {
		cout << endl << "RHS:" << endl;
		for(int i =0; i<nrow; i++) {
			cout << "( " << rows[i] << ", " << vals[i] << ")" <<endl;
		}

	}

	void finish_assembly() {}

	void apply_constrains(double) {}

	int solve() {}
};


TEST_F(SetValues, dirichlet) {
	vector<int> row_dofs { -1, -2, 3};
	vector<int> col_dofs { -1, -2, 3};
	arma::mat matrix { 0,2,3, 4,5,6, 7,8,9 };
	matrix.reshape(3,3);
	arma::vec rhs { 0,0,10 };
	arma::vec row_solution { 10, 20 , 100 };
	arma::vec col_solution { 10, 20 , 100 };

	this->set_values(row_dofs, col_dofs,matrix, rhs,
	        row_solution, col_solution);
}
