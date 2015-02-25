#define TEST_USE_MPI

#include <flow_gtest_mpi.hh>
#include "fields/vec_seq_double.hh"


TEST(VecSeqDouble, vec_data) {
	unsigned int data_size = 5;
	VectorSeqDouble v(data_size);
	typename VectorSeqDouble::VectorSeq data_ptr = v.get_data_ptr();

	EXPECT_EQ(data_size, data_ptr->size());
	EXPECT_DOUBLE_EQ(0.0, v[0]);
	v[0] = 2.5;
	EXPECT_DOUBLE_EQ(2.5, v[0]);
	EXPECT_DOUBLE_EQ((*data_ptr)[0], v[0]);
}
