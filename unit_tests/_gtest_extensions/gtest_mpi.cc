/*
 * gtest_flow.cc
 *
 *  Created on: Mar 15, 2012
 *      Author: jb
 */

#define TEST_USE_MPI        // include headers for MPI functions
#define TEST_HAS_MAIN       // do not include main in this unit
#include "gtest_mpi.hh"


namespace testing {
///////////////////////////////////////
// Following functions are unfortunately static in gtest.c, so we have no access to them BAD, BAD, BAD design

// Formats a countable noun.  Depending on its quantity, either the
// singular form or the plural form is used. e.g.
//
// FormatCountableNoun(1, "formula", "formuli") returns "1 formula".
// FormatCountableNoun(5, "book", "books") returns "5 books".
static std::string FormatCountableNoun(int count,
                                            const char * singular_form,
                                            const char * plural_form) {
	  return internal::StreamableToString(count) + " " +
	      (count == 1 ? singular_form : plural_form);
}

// Formats the count of tests.
static std::string FormatTestCount(int test_count) {
  return FormatCountableNoun(test_count, "test", "tests");
}

// Formats the count of test cases.
static std::string FormatTestCaseCount(int test_case_count) {
  return FormatCountableNoun(test_case_count, "test case", "test cases");
}

// Converts a TestPartResult::Type enum to human-friendly string
// representation.  Both kNonFatalFailure and kFatalFailure are translated
// to "Failure", as the user usually doesn't care about the difference
// between the two when viewing the test result.
static const char * TestPartResultTypeToString(TestPartResult::Type type) {
  switch (type) {
    case TestPartResult::kSuccess:
      return "Success";

    case TestPartResult::kNonFatalFailure:
    case TestPartResult::kFatalFailure:
#ifdef _MSC_VER
      return "error: ";
#else
      return "Failure\n";
#endif
    default:
      return "Unknown result type";
  }
}
// Prints a TestPartResult to a String.
static std::string PrintTestPartResultToString(
    const TestPartResult& test_part_result) {
  return (Message()
          << internal::FormatFileLocation(test_part_result.file_name(),
                                          test_part_result.line_number())
          << " " << TestPartResultTypeToString(test_part_result.type())
          << test_part_result.message()).GetString();
}

// Prints a TestPartResult.
static void PrintTestPartResult(const TestPartResult& test_part_result) {
  const internal::string& result =
      PrintTestPartResultToString(test_part_result);
  printf("%s\n", result.c_str());
  fflush(stdout);
  // If the test program runs in Visual Studio or a debugger, the
  // following statements add the test part result message to the Output
  // window such that the user can double-click on it to jump to the
  // corresponding source code location; otherwise they do nothing.
#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
  // We don't call OutputDebugString*() on Windows Mobile, as printing
  // to stdout is done by OutputDebugString() there already - we don't
  // want the same message printed twice.
  ::OutputDebugStringA(result.c_str());
  ::OutputDebugStringA("\n");
#endif
}


namespace internal {


MPI_PrettyUnitTestResultPrinter::MPI_PrettyUnitTestResultPrinter() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

}


void MPI_PrettyUnitTestResultPrinter::PrintTestName(const char * test_case, const char * test) {
  printf("%s.%s", test_case, test);
}



// Fired before each iteration of tests starts.
// (not yet MPI friendly)

void MPI_PrettyUnitTestResultPrinter::OnTestIterationStart(const UnitTest& unit_test, int iteration) {
}


void MPI_PrettyUnitTestResultPrinter::OnEnvironmentsSetUpStart(const UnitTest& /*unit_test*/) {
  if (rank==0) {
      ColoredPrintf(COLOR_GREEN,  "[----------] ");
      printf("Global test environment set-up.\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestCaseStart(const TestCase& test_case) {
  if (rank==0) {
      const std::string counts =
              FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");

      ColoredPrintf(COLOR_GREEN, "[----------] ");
      printf("%s from %s", counts.c_str(), test_case.name());
      if (test_case.type_param() == NULL) {
          printf("\n");
      } else {
          printf(", where TypeParam = %s\n", test_case.type_param());
      }
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestStart(const TestInfo& test_info) {
  if (rank == 0) {
      ColoredPrintf(COLOR_GREEN,  "[ RUN      ] ");
      PrintTestName(test_info.test_case_name(), test_info.name());
      printf("\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// Called after a failed assertion or a SUCCEED() invocation.
void MPI_PrettyUnitTestResultPrinter::OnTestPartResult(const TestPartResult& result) {
    // let print one process after other
    for(int i=0; i < np; i++) {
        if (rank == i && (result.type() != TestPartResult::kSuccess)) {
            printf("[%d] ", rank);
            PrintTestPartResult(result);
            fflush(stdout);
        }
    }
}

// (no yet fully MPI friendly)
void MPI_PrettyUnitTestResultPrinter::OnTestEnd(const TestInfo& test_info) {
  int success=(test_info.result()->Passed());
  int reduced_success;
  MPI_Allreduce(&success, &reduced_success, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);

  if (rank == 0) {
      if (reduced_success) {
          ColoredPrintf(COLOR_GREEN, "[       OK ] ");
      } else {
          ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
      }
      PrintTestName(test_info.test_case_name(), test_info.name());

      // HERE we should print comments for failed processes
      if (! reduced_success)
          PrintFullTestCommentIfPresent(test_info);

      if (GTEST_FLAG(print_time)) {
          printf(" (%s ms)\n", internal::StreamableToString(
                  test_info.result()->elapsed_time()).c_str());
      } else {
          printf("\n");
      }
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestCaseEnd(const TestCase& test_case) {
  if (!GTEST_FLAG(print_time)) return;

  if (rank==0) {
	  const std::string counts =
	      FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
	  ColoredPrintf(COLOR_GREEN, "[----------] ");
	  printf("%s from %s (%s ms total)\n\n",
	         counts.c_str(), test_case.name(),
	         internal::StreamableToString(test_case.elapsed_time()).c_str());
	  fflush(stdout);

  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnEnvironmentsTearDownStart(
    const UnitTest& /*unit_test*/) {
  if (rank==0) {
      ColoredPrintf(COLOR_GREEN,  "[----------] ");
      printf("Global test environment tear-down\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// Internal helper for printing the list of failed tests.
// (no yet MPI friendly)
void MPI_PrettyUnitTestResultPrinter::PrintFailedTests(const UnitTest& unit_test) {
  const int failed_test_count = unit_test.failed_test_count();
  if (failed_test_count == 0) {
    return;
  }

  for (int i = 0; i < unit_test.total_test_case_count(); ++i) {
    const TestCase& test_case = *unit_test.GetTestCase(i);
    if (!test_case.should_run() || (test_case.failed_test_count() == 0)) {
      continue;
    }
    for (int j = 0; j < test_case.total_test_count(); ++j) {
      const TestInfo& test_info = *test_case.GetTestInfo(j);
      if (!test_info.should_run() || test_info.result()->Passed()) {
        continue;
      }
      ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
      printf("%s.%s", test_case.name(), test_info.name());
      PrintFullTestCommentIfPresent(test_info);
      printf("\n");
    }
  }
}

// (no yet MPI friendly)
void MPI_PrettyUnitTestResultPrinter::OnTestIterationEnd(const UnitTest& unit_test,
                                                     int /*iteration*/) {
  ColoredPrintf(COLOR_GREEN,  "[==========] ");
  printf("%s from %s ran.",
         FormatTestCount(unit_test.test_to_run_count()).c_str(),
         FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
  if (GTEST_FLAG(print_time)) {
    printf(" (%s ms total)",
           internal::StreamableToString(unit_test.elapsed_time()).c_str());
  }
  printf("\n");
  ColoredPrintf(COLOR_GREEN,  "[  PASSED  ] ");
  printf("%s.\n", FormatTestCount(unit_test.successful_test_count()).c_str());

  int num_failures = unit_test.failed_test_count();
  if (!unit_test.Passed()) {
    const int failed_test_count = unit_test.failed_test_count();
    ColoredPrintf(COLOR_RED,  "[  FAILED  ] ");
    printf("%s, listed below:\n", FormatTestCount(failed_test_count).c_str());
    PrintFailedTests(unit_test);
    printf("\n%2d FAILED %s\n", num_failures,
                        num_failures == 1 ? "TEST" : "TESTS");
  }

  int num_disabled = unit_test.disabled_test_count();
  if (num_disabled && !GTEST_FLAG(also_run_disabled_tests)) {
    if (!num_failures) {
      printf("\n");  // Add a spacer if no FAILURE banner is displayed.
    }
    ColoredPrintf(COLOR_YELLOW,
                  "  YOU HAVE %d DISABLED %s\n\n",
                  num_disabled,
                  num_disabled == 1 ? "TEST" : "TESTS");
  }
  // Ensure that Google Test output is printed before, e.g., heapchecker output.
  fflush(stdout);
}

}} // end of testing and internal namespaces



