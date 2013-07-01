
#define DEBUG

#define TEST_USE_PETSC

#include <gtest_mpi.hh>

#include <la/schur.hh>

#include <petscmat.h>

/*TEST(la, matrix) {
	  Mat            A,A11,A12,A21,A22;
	  Vec            X,X1,X2,Y,Z,Z1,Z2;
	  PetscScalar    *a,*b,*x,*y,*z,v,one=1;
	  PetscReal      nrm;
	  PetscErrorCode ierr;
	  PetscInt       size=8,size1=6,size2=2, i,j;

	  //PetscInitialize(&argc,&argv,0,help);

	  // Create matrix and three vectors: these are all normal
	  ierr = PetscMalloc(size*size*sizeof(PetscScalar),&a);
	  ierr = PetscMalloc(size*size*sizeof(PetscScalar),&b);
	  for (i=0; i<size; i++) {
	    for (j=0; j<size; j++) {
	      a[i+j*size] = rand(); b[i+j*size] = a[i+j*size];
	    }
	  }
	  ierr = MatCreate(MPI_COMM_SELF,&A);
	  ierr = MatSetSizes(A,size,size,size,size);
	  ierr = MatSetType(A,MATSEQDENSE);
	  ierr = MatSeqDenseSetPreallocation(A,a);
	  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	  ierr = PetscMalloc(size*sizeof(PetscScalar),&x);
	  for (i=0; i<size; i++) {
	    x[i] = one;
	  }
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size,x,&X);
	  ierr = VecAssemblyBegin(X);
	  ierr = VecAssemblyEnd(X);

	  ierr = PetscMalloc(size*sizeof(PetscScalar),&y);
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size,y,&Y);
	  ierr = VecAssemblyBegin(Y);
	  ierr = VecAssemblyEnd(Y);

	  ierr = PetscMalloc(size*sizeof(PetscScalar),&z);
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size,z,&Z);
	  ierr = VecAssemblyBegin(Z);
	  ierr = VecAssemblyEnd(Z);

	  // Now create submatrices and subvectors
	  ierr = MatCreate(MPI_COMM_SELF,&A11);
	  ierr = MatSetSizes(A11,size1,size1,size1,size1);
	  ierr = MatSetType(A11,MATSEQDENSE);
	  ierr = MatSeqDenseSetPreallocation(A11,b);
	  ierr = MatSeqDenseSetLDA(A11,size);
	  ierr = MatAssemblyBegin(A11,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A11,MAT_FINAL_ASSEMBLY);

	  ierr = MatCreate(MPI_COMM_SELF,&A12);
	  ierr = MatSetSizes(A12,size1,size2,size1,size2);
	  ierr = MatSetType(A12,MATSEQDENSE);
	  ierr = MatSeqDenseSetPreallocation(A12,b+size1*size);
	  ierr = MatSeqDenseSetLDA(A12,size);
	  ierr = MatAssemblyBegin(A12,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A12,MAT_FINAL_ASSEMBLY);

	  ierr = MatCreate(MPI_COMM_SELF,&A21);
	  ierr = MatSetSizes(A21,size2,size1,size2,size1);
	  ierr = MatSetType(A21,MATSEQDENSE);
	  ierr = MatSeqDenseSetPreallocation(A21,b+size1);
	  ierr = MatSeqDenseSetLDA(A21,size);
	  ierr = MatAssemblyBegin(A21,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A21,MAT_FINAL_ASSEMBLY);

	  ierr = MatCreate(MPI_COMM_SELF,&A22);
	  ierr = MatSetSizes(A22,size2,size2,size2,size2);
	  ierr = MatSetType(A22,MATSEQDENSE);
	  ierr = MatSeqDenseSetPreallocation(A22,b+size1*size+size1);
	  ierr = MatSeqDenseSetLDA(A22,size);
	  ierr = MatAssemblyBegin(A22,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A22,MAT_FINAL_ASSEMBLY);

	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size1,x,&X1);
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size2,x+size1,&X2);
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size1,z,&Z1);
	  ierr = VecCreateSeqWithArray(MPI_COMM_SELF,1,size2,z+size1,&Z2);

	  // Now multiple matrix times input in two ways;
	  // compare the results
	  ierr = MatMult(A,X,Y);
	  ierr = MatMult(A11,X1,Z1);
	  ierr = MatMultAdd(A12,X2,Z1,Z1);
	  ierr = MatMult(A22,X2,Z2);
	  ierr = MatMultAdd(A21,X1,Z2,Z2);
	  ierr = VecAXPY(Z,-1.0,Y);
	  ierr = VecNorm(Z,NORM_2,&nrm);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test1; error norm=%G\n",nrm);

	  ierr = PetscPrintf(PETSC_COMM_WORLD,"MatMult the usual way:\n");
	  ierr = VecView(Y,0);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"MatMult by subblock:\n");
	  ierr = VecView(Z,0);

	  // Next test: change both matrices
	  v    = rand(); i=1; j=size-2;
	  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);
	  j   -= size1;
	  ierr = MatSetValues(A12,1,&i,1,&j,&v,INSERT_VALUES);
	  v    = rand(); i=j=size1+1;
	  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);
	  i    =j=1;
	  ierr = MatSetValues(A22,1,&i,1,&j,&v,INSERT_VALUES);
	  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyBegin(A12,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A12,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyBegin(A22,MAT_FINAL_ASSEMBLY);
	  ierr = MatAssemblyEnd(A22,MAT_FINAL_ASSEMBLY);

	  ierr = MatMult(A,X,Y);
	  ierr = MatMult(A11,X1,Z1);
	  ierr = MatMultAdd(A12,X2,Z1,Z1);
	  ierr = MatMult(A22,X2,Z2);
	  ierr = MatMultAdd(A21,X1,Z2,Z2);
	  ierr = VecAXPY(Z,-1.0,Y);
	  ierr = VecNorm(Z,NORM_2,&nrm);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test2; error norm=%G\n",nrm);

	  // Transpose product
	  ierr = MatMultTranspose(A,X,Y);
	  ierr = MatMultTranspose(A11,X1,Z1);
	  ierr = MatMultTransposeAdd(A21,X2,Z1,Z1);
	  ierr = MatMultTranspose(A22,X2,Z2);
	  ierr = MatMultTransposeAdd(A12,X1,Z2,Z2);
	  ierr = VecAXPY(Z,-1.0,Y);
	  ierr = VecNorm(Z,NORM_2,&nrm);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Test3; error norm=%G\n",nrm);

	  ierr = PetscFree(a);
	  ierr = PetscFree(b);
	  ierr = PetscFree(x);
	  ierr = PetscFree(y);
	  ierr = PetscFree(z);
	  ierr = MatDestroy(&A);
	  ierr = MatDestroy(&A11);
	  ierr = MatDestroy(&A12);
	  ierr = MatDestroy(&A21);
	  ierr = MatDestroy(&A22);

	  ierr = VecDestroy(&X);
	  ierr = VecDestroy(&Y);
	  ierr = VecDestroy(&Z);

	  ierr = VecDestroy(&X1);
	  ierr = VecDestroy(&X2);
	  ierr = VecDestroy(&Z1);
	  ierr = VecDestroy(&Z2);

	  ierr = PetscFinalize();
} // */

TEST(la, inversion_matrix) {
	Mat matrix;
	PetscInt size = 15;
	PetscErrorCode ierr;
	PetscScalar *a;

	ierr = PetscMalloc(size*size*sizeof(PetscScalar),&a);
	for (PetscInt i=0; i<size; i++) {
		for (PetscInt j=0; j<size; j++) {
			a[i+j*size] = 0;
		}
    }

	ierr = MatCreate(MPI_COMM_SELF,&matrix);
	ierr = MatSetSizes(matrix,size,size,size,size);
	ierr = MatSetType(matrix,MATAIJ);
	ierr = MatSeqDenseSetPreallocation(matrix, a);
	ierr = MatSetUp(matrix);

	PetscInt row [] =   {0,0,0, 0,1,1,1,2,2, 2,3,3,3,4,4,5,5,5,6,6,7,8,8,9,9,9,10,10,10,11,11,12,13,13,13,14,14};
	PetscInt col [] =   {1,2,3,10,0,2,3,1,2,10,0,2,3,4,5,4,8,9,6,7,6,5,8,4,8,9, 0, 1,10,13,14,13,11,12,14,11,12};
	PetscScalar val[] = {1,2,3, 4,5,6,7,8,9,10,1,2,3,1,2,3,4,5,3,2,1,5,6,7,8,9,11,12,13, 1, 2, 3, 4, 5, 6, 7, 8};

	for (int i = 0; i < 37; i++) {
		ierr = MatSetValue(matrix, row[i], col[i], val[i], INSERT_VALUES);
	}

	// test output
	/*PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;
	for (PetscInt i=0; i<size; i++) {
		ierr = MatGetRow(matrix, i, &ncols, &cols, &vals);
		printf("line %d: ", i );
		for (PetscInt j=0; j<size; j++) {
			printf("%d->%f ", (int) (cols[j]), (double) (vals[j]) );
		}
		printf("\n");
	} // */

	ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);

	SchurComplement schurComplement(matrix, 20);

}
