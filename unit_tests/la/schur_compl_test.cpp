

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

class SchurComplementTest : public SchurComplement {
public:
	SchurComplementTest(IS ia, Distribution *ds)
	: SchurComplement(ia, ds)
	{}

	Mat get_a_inv() const {return (IA);}

	/**
	 * Fill random local part of block matrix
	 * A  B
	 * Bt 0
	 *
	 * where A is block diagonal. Local blocks sizes are rows[min_idx] .. rows[max_idx-1].
	 * Block B has number of columns equal to number of blocks.
	 */
	void fill_matrix(int rank, Distribution &ds, Distribution &block_ds) {

		// set B columns
		int n_cols_B=block_ds.size();
		std::vector<PetscInt> b_cols(n_cols_B);
		for( unsigned int p=0;p<block_ds.np();p++)
			for (unsigned int j=block_ds.begin(p); j<block_ds.end(p); j++) {
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
			mat_set_values(block_size, &a_rows[0], block_size, &a_rows[0], &a_vals[0]);

			// set B values
			std::vector<PetscScalar> b_vals(block_size*n_cols_B);
			for (int j=0; j<block_size*n_cols_B; j++)
				b_vals[j] = 1;

			// set C values
			std::vector<PetscScalar> c_vals(n_cols_B);
			for (int j=0; j<n_cols_B; j++)
				c_vals[j] = 0;

			// must iterate per rows to get correct transpose
			for(unsigned int row=0; row<block_size;row++) {
				mat_set_values(1, &a_rows[row], 1, &b_cols[rank], &b_vals[row*n_cols_B]);
				mat_set_values(1, &b_cols[rank],1, &a_rows[row], &b_vals[row*n_cols_B]);
			}

			mat_set_values(1, &b_cols[rank], 1, &b_cols[rank], &c_vals[rank]);

		}
	}
};

class LinSysPetscTest : public LinSys_PETSC {
public:
	LinSysPetscTest(Distribution *ds)
	: LinSys_PETSC(ds)
	{ r_tol_ = 1e-12; a_tol_ = 1e-12; }
};


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
	SchurComplementTest * schurComplement = new SchurComplementTest(set, &all_ds);
	schurComplement->set_solution(NULL);
	schurComplement->set_positive_definite();
	schurComplement->start_allocation();
	schurComplement->fill_matrix( rank, ds, block_ds); // preallocate matrix
	schurComplement->start_add_assembly();
	schurComplement->fill_matrix( rank, ds, block_ds); // fill matrix
	schurComplement->finish_assembly();
	MatView(*(schurComplement->get_matrix()),PETSC_VIEWER_STDOUT_WORLD);

	LinSys * lin_sys = new LinSysPetscTest( schurComplement->make_complement_distribution() );
	schurComplement->set_complement( (LinSys_PETSC *)lin_sys );
	schurComplement->solve();

	// test of computed values
	{
		PetscInt ncols;
		const PetscInt *cols;
		const PetscScalar *vals;
		for (unsigned int i=0; i<block_size; i++) {
			MatGetRow(schurComplement->get_a_inv(), i + rank*block_size, &ncols, &cols, &vals);
			EXPECT_FLOAT_EQ( (1.0 / (double)(rank + 2)), vals[i] );
			MatRestoreRow(schurComplement->get_a_inv(), i + rank*block_size, &ncols, &cols, &vals);
		}
		MatGetRow(*(schurComplement->get_system()->get_matrix()), rank, &ncols, &cols, &vals);
		EXPECT_FLOAT_EQ( ((double)block_size / (double)(rank + 2)), vals[0] );
		MatRestoreRow(*(schurComplement->get_system()->get_matrix()), rank, &ncols, &cols, &vals);
	}
}
