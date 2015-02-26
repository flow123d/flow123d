#define TEST_USE_MPI

#include <flow_gtest_mpi.hh>
#include "fields/vec_seq_double.hh"


TEST(VecSeqDouble, vec_data) {
	PetscInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

	unsigned int data_size = 5;
	VectorSeqDouble v;
	v.resize(data_size);
	typename VectorSeqDouble::VectorSeq data_ptr = v.get_data_ptr();

	EXPECT_EQ(data_size, data_ptr->size());
	EXPECT_DOUBLE_EQ(0.0, v[0]);
	v[0] = 2.5;
	EXPECT_DOUBLE_EQ(2.5, v[0]);
	EXPECT_DOUBLE_EQ((*data_ptr)[0], v[0]);

	Vec petscVec = v.get_data_petsc();
	unsigned int indices[5] = { 0,1,2,3,4 };
	double vals[5];
	VecGetValues(petscVec, data_size, (PetscInt *)indices, vals);
	EXPECT_DOUBLE_EQ(2.5, vals[0]);
	EXPECT_DOUBLE_EQ(0.0, vals[1]);
}
