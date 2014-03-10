#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

TEST(trivial, trivial) {
  EXPECT_EQ(1,1);
}