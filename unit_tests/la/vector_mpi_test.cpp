#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>
#include "la/vector_mpi.hh"


TEST(VecMPI, vec_data) {
	PetscInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    int nproc;
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    if (nproc != 1) return; // continue only in sequential mode

	unsigned int data_size = 5;
	VectorMPI v(data_size);

	EXPECT_EQ(data_size, v.size());
	EXPECT_DOUBLE_EQ(0.0, v.get(0));
	v[0] = 2.5;
	EXPECT_DOUBLE_EQ(2.5, v.get(0));
	EXPECT_DOUBLE_EQ(v.get(0), v[0]);

	Vec petscVec = v.petsc_vec();
	unsigned int indices[5] = { 0,1,2,3,4 };
	double vals[5];
	VecGetValues(petscVec, data_size, (PetscInt *)indices, vals);
	EXPECT_DOUBLE_EQ(2.5, vals[0]);
	EXPECT_DOUBLE_EQ(0.0, vals[1]);

	v.zero_entries();
	EXPECT_DOUBLE_EQ(0.0, v[0]);
}



TEST(VecMPI, ghost_values) {
    PetscInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    
    int rank, nproc;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    if (nproc == 1) return; // continue only in parallel mode
    
    unsigned int data_size = 1;
    vector<LongIdx> ghost_idx = { ((rank+1)%nproc) };
    VectorMPI v(data_size, ghost_idx);
    EXPECT_EQ( data_size+ghost_idx.size(), v.size() );
    
    v[0] = 1;
    v.local_to_ghost_begin();
    v.local_to_ghost_end();
    EXPECT_DOUBLE_EQ( 1, v[data_size] ); // ghost value should be equal to the local value from proc. rank+1
    
    v[data_size] = 2;
    v.ghost_to_local_begin();
    v.ghost_to_local_end();
    EXPECT_DOUBLE_EQ( 3, v[0] ); // ghost values are added to the local value
}
