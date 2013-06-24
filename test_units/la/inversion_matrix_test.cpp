
#define DEBUG

#define TEST_USE_PETSC

#include <gtest_mpi.hh>

#include <la/schur.hh>

//#include "petscmat.h"

TEST(inversion, inversion_matrix) {
	Mat matrix;
	PetscInt size = 15;

	MatCreate(MPI_COMM_SELF, &matrix);
	MatSetSizes(matrix, size, size, size, size);

	PetscInt row [] =   {0,0,0, 0,1,1,1,2,2, 2,3,3,3,4,4,5,5,5,6,6,7,8,8,9,9,9,10,10,10,11,11,12,13,13,13,14,14};
	PetscInt col [] =   {1,2,3,10,0,2,3,1,2,10,0,2,3,4,5,4,8,9,6,7,6,5,8,4,8,9, 0, 1,10,13,14,13,11,12,14,11,12};
	PetscScalar val[] = {1,2,3, 4,5,6,7,8,9,10,1,2,3,1,2,3,4,5,3,2,1,5,6,7,8,9,11,12,13, 1, 2, 3, 4, 5, 6, 7, 8};

	for (int i = 0; i < 37; i++) MatSetValue(matrix, row[i], col[i], val[i], INSERT_VALUES);

	SchurComplement schurComplement(matrix);

}
