

#define TEST_USE_PETSC

#include "flow_gtest_mpi.hh"

#include <la/distribution.hh>
#include <la/schur.hh>
#include <la/linsys.hh>
#include "la/linsys_PETSC.hh"

#include <petscmat.h>
#include <math.h>
#include <vector>


const int block_size = 2;
const int block_count = 1;

/**
 * Fill random local part of block matrix
 * A  B
 * Bt 0
 *
 * where A is block diagonal. Local blocks sizes are rows[min_idx] .. rows[max_idx-1].
 * Block B has number of columns equal to number of blocks.
 */
/*void fill_matrix(LinSys * lin_sys, int *blocks, Distribution &ds, Distribution &block_ds) {


	// set B columns
	int n_cols_B=block_ds.size();
	std::vector<PetscInt> b_cols(n_cols_B);
	for( int p=0;p<block_ds.np();p++)
		for (unsigned int j=block_ds.begin(p); j<block_ds.end(p); j++) {
			int proc=block_ds.get_proc(j);
			b_cols[j]=ds.end(p)+j;
		}


	// create block A of matrix
	int local_idx=0;
	for (unsigned int i = block_ds.begin(); i < block_ds.end(); i++) {
		int block_size=blocks[i];
		// make random block values
		std::vector<PetscScalar> a_vals(block_size * block_size);
		for (unsigned int j=0; j<block_size*block_size; j++)
			a_vals[j]= rand()%19 + 1;

		// set rows and columns indices
		std::vector<PetscInt> a_rows(block_size);
		for (unsigned int j=0; j<block_size; j++) {
			a_rows[j]=ds.begin() + block_ds.begin() + local_idx;
			local_idx++;
		}
		lin_sys->mat_set_values(block_size, &a_rows[0], block_size, &a_rows[0], &a_vals[0]);

		// set B values
		std::vector<PetscScalar> b_vals(block_size*n_cols_B);
		for (unsigned int j=0; j<block_size*n_cols_B; j++)
			b_vals[j] =rand()%19 + 1;

		// must iterate per rows to get correct transpose
		for(unsigned int row=0; row<block_size;row++) {
			lin_sys->mat_set_values(1, &a_rows[row], n_cols_B, &b_cols[0], &b_vals[row*n_cols_B]);
			lin_sys->mat_set_values(n_cols_B, &b_cols[0],1, &a_rows[row], &b_vals[row*n_cols_B]);
		}

	}
}*/

void fill_matrix(LinSys * lin_sys, int rank, Distribution &ds, Distribution &block_ds) {

	// set B columns
	int n_cols_B=block_ds.size();
	std::vector<PetscInt> b_cols(n_cols_B);
	for( int p=0;p<block_ds.np();p++)
		for (unsigned int j=block_ds.begin(p); j<block_ds.end(p); j++) {
			int proc=block_ds.get_proc(j);
			b_cols[j]=ds.end(p)+j;
		}

	// create block A of matrix
	int local_idx=0;
	for (unsigned int i = block_ds.begin(); i < block_ds.end(); i++) {
		// make random block values
		std::vector<PetscScalar> a_vals(block_size * block_size, 0);
		for (unsigned int j=0; j<block_size; j++)
			a_vals[ j + j*block_size ]= (rank + 2);

		// set rows and columns indices
		std::vector<PetscInt> a_rows(block_size);
		for (unsigned int j=0; j<block_size; j++) {
			a_rows[j]=ds.begin() + block_ds.begin() + local_idx;
			local_idx++;
		}
		lin_sys->mat_set_values(block_size, &a_rows[0], block_size, &a_rows[0], &a_vals[0]);

		// set B values
		std::vector<PetscScalar> b_vals(block_size*n_cols_B);
		for (unsigned int j=0; j<block_size*n_cols_B; j++)
			b_vals[j] = 1;

		// must iterate per rows to get correct transpose
		for(unsigned int row=0; row<block_size;row++) {
			lin_sys->mat_set_values(1, &a_rows[row], n_cols_B, &b_cols[0], &b_vals[row*n_cols_B]);
			lin_sys->mat_set_values(n_cols_B, &b_cols[0],1, &a_rows[row], &b_vals[row*n_cols_B]);
		}

	}
}

TEST(schur, complement) {
	IS set;
	// vytvorit rozdeleni bloku na procesory ve tvaru "part" (tj. indexy prvnich radku na procesorech)
    int np, rank;

    MPI_Comm_size(PETSC_COMM_WORLD, &np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    Distribution ds(block_size, MPI_COMM_WORLD);
    Distribution block_ds(block_count, MPI_COMM_WORLD);
    Distribution all_ds(block_size + block_count, MPI_COMM_WORLD);
    /*if (rank == 0) {
        cout << all_ds;
        cout << ds;
        cout << block_ds;
    }*/

	ISCreateStride(PETSC_COMM_WORLD, ds.lsize(), all_ds.begin(), 1, &set);
	ISView(set, PETSC_VIEWER_STDOUT_WORLD);

    // volat s lokalni velkosti = pocet radku na lokalnim proc.
	SchurComplement * schurComplement = new SchurComplement(set, &all_ds);
	schurComplement->set_solution(NULL);
	schurComplement->set_symmetric();
	schurComplement->start_allocation();
	fill_matrix( schurComplement, rank, ds, block_ds); // preallocate matrix
	VecZeroEntries(schurComplement->get_solution());
	schurComplement->start_add_assembly();
	VecZeroEntries(schurComplement->get_solution());
	fill_matrix( schurComplement, rank, ds, block_ds); // fill matrix
	schurComplement->finish_assembly();
	MatView(schurComplement->get_matrix(),PETSC_VIEWER_STDOUT_WORLD);

	LinSys * lin_sys = new LinSys_PETSC( schurComplement->make_complement_distribution() );
	schurComplement->set_complement( (LinSys_PETSC *)lin_sys );
	schurComplement->create_inversion_matrix();
	MatView(schurComplement->get_a_inv(),PETSC_VIEWER_STDOUT_WORLD);
	schurComplement->form_schur();
	//schurComplement->set_spd();

}

/*TEST(schur, inversion_matrix) {
	int blocks [] = {5,2,3,3,4,2};
	int n_blocks = 6;
	int max_block_size=5;



	int first_idx=0, size=0;
	IS set;
	// vytvorit rozdeleni bloku na procesory ve tvaru "part" (tj. indexy prvnich radku na procesorech)
    int np, rank;
    double block_size;
    int min_idx, max_idx;

    MPI_Comm_size(PETSC_COMM_WORLD, &np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// assign blocks to processors
    int total_size=0;
    for(int i=0;i<n_blocks;i++) total_size+=blocks[i];
    vector<int> blocks_part(n_blocks);
    int proc=0;
    int distributed=0;
    int local_size=0;
    int block_local_size=0;
    for(int i=0;i<n_blocks;i++) {
    	blocks_part[i]=proc;
    	distributed+=blocks[i];
    	if (proc==rank) {
    		local_size+=blocks[i];
    		block_local_size++;
    	}
    	if (distributed > total_size*(proc+1)/np) proc++;

    }

    Distribution ds(local_size, MPI_COMM_WORLD);
    Distribution block_ds(block_local_size, MPI_COMM_WORLD);
    Distribution all_ds(local_size+block_local_size, MPI_COMM_WORLD);
    //cout << ds;
    //cout << block_ds;

    // volat s lokalni velkosti = pocet radku na lokalnim proc.
	LinSys * lin_sys = new LinSys_PETSC(&all_ds);
	lin_sys->set_solution(NULL);
	lin_sys->set_symmetric();
	lin_sys->start_allocation();
	time_t seed=time(NULL);
	srand(seed);
	fill_matrix( lin_sys, blocks, ds, block_ds); // preallocate matrix
	lin_sys->start_add_assembly();
	srand(seed);
	fill_matrix( lin_sys, blocks, ds, block_ds); // fill matrix
	lin_sys->finish_assembly();
	MatView(lin_sys->get_matrix(),PETSC_VIEWER_STDOUT_WORLD);

	ISCreateStride(PETSC_COMM_WORLD, ds.lsize(), all_ds.begin(), 1, &set);
	ISView(set, PETSC_VIEWER_STDOUT_WORLD);

	SchurComplement schurComplement(lin_sys, set, &ds);

}*/
