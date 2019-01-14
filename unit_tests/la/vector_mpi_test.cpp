#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>
#include "la/vector_mpi.hh"


TEST(VecMPI, vec_data) {
	PetscInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

	unsigned int data_size = 5;
	VectorMPI v(data_size);
	typename VectorMPI::VectorDataPtr data_ptr = v.data_ptr();

	EXPECT_EQ(data_size, data_ptr->size());
	EXPECT_DOUBLE_EQ(0.0, v[0]);
	v[0] = 2.5;
	EXPECT_DOUBLE_EQ(2.5, v[0]);
	EXPECT_DOUBLE_EQ((*data_ptr)[0], v[0]);

	Vec petscVec = v.petsc_vec();
	unsigned int indices[5] = { 0,1,2,3,4 };
	double vals[5];
	VecGetValues(petscVec, data_size, (PetscInt *)indices, vals);
	EXPECT_DOUBLE_EQ(2.5, vals[0]);
	EXPECT_DOUBLE_EQ(0.0, vals[1]);

	v.zero_entries();
	EXPECT_DOUBLE_EQ(0.0, v[0]);
}
