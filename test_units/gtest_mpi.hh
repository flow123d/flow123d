/*
 * gtest_flow.hh
 *
 *  Created on: Mar 14, 2012
 *      Author: jb
 */

/**
 * For tests of flow123d it should be enough to include this file and files
 * that will be tested. Before include you can define some macros to influence
 * functionality provided by gtest_flow.hh.
 *
 * TEST_HAS_MAIN - do not provide main function, if not specified gtest_flow.hh fabricates its own main dunction compatible
 *                 with possible usage of MPI and PETSc
 *
 * TEST_USE_PETSc - Initialize PETSc and include petsc.h. You has to include other modules yourself. Set automatically
 *                  TEST_USE_MPI.
 *
 * TEST_USE_MPI -   Implements Listener with support of parallel applications. Output of individual processes
 *                  is sorted according to rank within one test. Result of a test is reduced over all processors, so the
 *                  test succeed only if all processes succeed. Header of test is reported only once.
 */

#ifndef GTEST_FLOW_HH_
#define GTEST_FLOW_HH_

//#include <global_defs.h>
//#include <system/system.hh>

#include <gtest/gtest.h>

#ifdef TEST_USE_PETSC
#define TEST_USE_MPI
#include <petsc.h>
#endif

#ifdef TEST_USE_MPI
#include <mpi.h>
#endif


#ifdef TEST_USE_MPI

/*
 * Implement modification of default event listener PrettyUnitTestResultPrinter that is
 * MPI compatible.
 */
namespace testing {

namespace internal {

//#include "gtest/internal/gtest-port.h"

// various imported thinks
// A test filter that matches everything.
static const char kUniversalFilter[] = "*";

enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

void ColoredPrintf(GTestColor color, const char* fmt, ...);
void PrintFullTestCommentIfPresent(const TestInfo& test_info);

// This class implements the TestEventListener interface.
//
// Class PrettyUnitTestResultPrinter is copyable.
class MPI_PrettyUnitTestResultPrinter : public ::testing::TestEventListener {
public:
  // Get basic MPI info.
  MPI_PrettyUnitTestResultPrinter() {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &np);

  }

  static void PrintTestName(const char * test_case, const char * test) {
    printf("%s.%s", test_case, test);
  }

  // The following methods override what's in the TestEventListener class.
  virtual void OnTestProgramStart(const UnitTest& /*unit_test*/) {}
  virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
  virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsSetUpEnd(const UnitTest& /*unit_test*/) {}
  virtual void OnTestCaseStart(const TestCase& test_case);
  virtual void OnTestStart(const TestInfo& test_info);
  virtual void OnTestPartResult(const TestPartResult& result);
  virtual void OnTestEnd(const TestInfo& test_info);
  virtual void OnTestCaseEnd(const TestCase& test_case);
  virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsTearDownEnd(const UnitTest& /*unit_test*/) {}
  virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
  virtual void OnTestProgramEnd(const UnitTest& /*unit_test*/) {}

private:
  int np, rank;


  static void PrintFailedTests(const UnitTest& unit_test);
  internal::String test_case_name_;
};


} // internal namespace
} // testing namespace
#endif // TEST_USE_MPI

#ifndef TEST_HAS_MAIN
int main(int argc, char** argv) {


  // This allows the user to override the flag on the command line.
  ::testing::InitGoogleTest(&argc, argv);
#ifdef TEST_USE_MPI
#ifdef TEST_USE_PETSC
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
#else
  MPI_Init(&argc, &argv);
#endif


  // Gets hold of the event listener list.
  ::testing::TestEventListeners& listeners =  ::testing::UnitTest::GetInstance()->listeners();
  delete listeners.Release(listeners.default_result_printer());

  // Adds a listener to the end.  Google Test takes the ownership.
  listeners.Append(new ::testing::internal::MPI_PrettyUnitTestResultPrinter);
#endif

  return RUN_ALL_TESTS();

#ifdef TEST_USE_MPI
#ifdef TEST_USE_PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif
#endif
}

#endif // not def TEST_HAS_MAIN




#endif /* GTEST_FLOW_HH_ */
