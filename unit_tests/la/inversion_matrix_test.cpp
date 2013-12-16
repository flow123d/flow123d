

#define TEST_USE_PETSC

#include <flow_gtest_mpi.hh>

#include <la/schur.hh>
#include <la/linsys.hh>

#include <petscmat.h>
#include <math.h>


PetscInt rows [] = {5,2,3,3,4,2};


void fill_matrix(LinSys * lin_sys, int min_idx, int max_idx) {
	unsigned int row_cnt = 0; //row counter
	unsigned int local_min, local_max;
	for (unsigned int i = 0; i < min_idx; i++) {
		row_cnt += rows[i];
	}
	local_min = row_cnt;

	// create block A of matrix
	for (unsigned int i = min_idx; i < max_idx; i++) {
		std::vector<PetscScalar> a_vals;
		a_vals.reserve(rows[i] * rows[i]);
		for (unsigned int j=0; j<rows[i]*rows[i]; j++) {
			a_vals.push_back( rand()%19 + 1 );
		}
		std::vector<int> a_rows_idx;
		a_rows_idx.reserve(rows[i]);
		for (unsigned int j=0; j<rows[i]; j++) {
			a_rows_idx.push_back(row_cnt + min_idx + j);
		}
		int * insert_rows = &a_rows_idx[0];
		PetscScalar * insert_vals = &a_vals[0];
		lin_sys->mat_set_values(rows[i], insert_rows, rows[i], insert_rows, insert_vals);
		row_cnt += rows[i];
	}
	local_max = row_cnt;
	/*for (unsigned int i = max_idx; i < 6; i++) {
		row_cnt += rows[i];
	}*/

	// create block B of matrix
	std::vector<PetscScalar> b_vals((local_max - local_min) * (max_idx - min_idx), 1.0);
	std::vector<int> b_rows_idx;
	std::vector<int> b_cols_idx;
	b_rows_idx.reserve(local_max - local_min);
	b_cols_idx.reserve(max_idx - min_idx);
	for (unsigned int i = local_min; i < local_max; i++) {
		b_rows_idx.push_back(i + min_idx);
	}
	for (unsigned int i = min_idx; i < max_idx; i++) {
		b_cols_idx.push_back(local_max + i);
	}
	int * b_insert_rows = &b_rows_idx[0];
	int * b_insert_cols = &b_cols_idx[0];
	PetscScalar * b_insert_vals = &b_vals[0];
	lin_sys->mat_set_values((local_max - local_min), b_insert_rows, (max_idx - min_idx), b_insert_cols, b_insert_vals);
	lin_sys->mat_set_values((max_idx - min_idx), b_insert_cols, (local_max - local_min), b_insert_rows, b_insert_vals); // */
}


TEST(la, inversion_matrix) {
	srand(123);

	int first_idx=0, size=0, submat_blocks=6;
	IS set;
	// vytvorit rozdeleni bloku na procesory ve tvaru "part" (tj. indexy prvnich radku na procesorech)
    int np, rank;
    double block_size;
    int min_idx, max_idx;

    MPI_Comm_size(PETSC_COMM_WORLD, &np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    block_size = (double)submat_blocks / (double)np;
    min_idx = (int) ( round(block_size * rank) );
    max_idx = (int) ( round(block_size * (rank + 1)) );
    for (int i = 0; i < min_idx; i++) {
    	first_idx += rows[i];
    }
    for (int i = min_idx; i < max_idx; i++) {
    	size += rows[i];
    }

    // volat s lokalni velkosti = pocet radku na lokalnim proc.
	LinSys * lin_sys = new LinSys_MPIAIJ(size + max_idx - min_idx);
	lin_sys->set_symmetric();
	lin_sys->start_allocation();
	fill_matrix( lin_sys, min_idx, max_idx ); // preallocate matrix
	lin_sys->start_add_assembly();
	fill_matrix( lin_sys, min_idx, max_idx ); // fill matrix
	lin_sys->finalize();
	MatView(lin_sys->get_matrix(),PETSC_VIEWER_STDOUT_WORLD);

	ISCreateStride(PETSC_COMM_WORLD, size, first_idx + min_idx, 1, &set); // kazdy proc. lokalni cast indexsetu viz. schur.cc line 386
	ISView(set, PETSC_VIEWER_STDOUT_WORLD);

	SchurComplement schurComplement(lin_sys, set, 6);

}
