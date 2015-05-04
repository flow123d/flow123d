#define TEST_USE_PETSC
#define TEST_HAS_MAIN

#include <flow_gtest_mpi.hh>
#include "deal.II/include/deal.II/lac/vector.templates.h"
#include "deal.II/include/deal.II/lac/sparse_direct.h"
#include "la/distribution.hh"


TEST(vector, sequential)
{
	// Try to create and write to a sequential vector.

	dealii::Vector<double> vec_seq(10);
	vec_seq = 1;

	EXPECT_DOUBLE_EQ( 10, vec_seq.norm_sqr() );

}

TEST(vector, petsc)
{
	// Try to create PETSc vector.

	dealii::PETScWrappers::Vector vec_petsc(10);
	vec_petsc(0) = 1;
}

TEST(vector, petsc_mpi)
{
	// Try to create PETSc parallel vector.

	Distribution ds(DistributionBlock(), 100, PETSC_COMM_WORLD);
	dealii::PETScWrappers::MPI::Vector vec_petsc(ds.get_comm(), ds.size(), ds.lsize());
	vec_petsc(0) = 1;
}

TEST(solver, sparse_direct_umfpack)
{
	// Test UMFPACK solver on the system
	//
	// ( 1 1 ) ( x ) = ( 5 )
	// ( 1 0 ) ( y ) = ( 2 )
	//
	// with the solution x = 2, y = 3.

	dealii::SparseDirectUMFPACK solver;
	dealii::SparsityPattern pattern(2, 2);
	dealii::SparseMatrix<double> mat;
	dealii::Vector<double> rhs(2);

	pattern.add(0, 1);
	pattern.add(1, 0);
	pattern.compress();
	mat.reinit(pattern);

	mat.set(0, 0, 1);
	mat.set(0, 1, 1);
	mat.set(1, 0, 1);

	rhs(0) = 5;
	rhs(1) = 2;
	
	solver.solve(mat, rhs);
	
	EXPECT_DOUBLE_EQ( 2, rhs(0) );
	EXPECT_DOUBLE_EQ( 3, rhs(1) );
}




int main(int argc, char** argv)
{
	// This is necessary when deal.II is built with threads.
	dealii::multithread_info.set_thread_limit(1);

	// This allows the user to override the flag on the command line.
	::testing::InitGoogleTest(&argc, argv);
	PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);


	// Gets hold of the event listener list.
	::testing::TestEventListeners& listeners =  ::testing::UnitTest::GetInstance()->listeners();
	delete listeners.Release(listeners.default_result_printer());

	// Adds a listener to the end.  Google Test takes the ownership.
	listeners.Append(new ::testing::internal::MPI_PrettyUnitTestResultPrinter);

	return RUN_ALL_TESTS();

	PetscFinalize();
}
