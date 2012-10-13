/*
 * gtest_throw_what.hh
 *
 *  Created on: May 21, 2012
 *      Author: jb
 */

#ifndef GTEST_THROW_WHAT_HH_
#define GTEST_THROW_WHAT_HH_

/**
 * Macros for death test that test particular class of exception
 * and more over test contents of message created by its what() method.
 */


#include <gtest/gtest.h>

// Returns an indented copy of stderr output for a death test.
// This makes distinguishing death test output lines from regular log lines
// much easier.
static ::std::string FormatDeathTestOutput(const ::std::string& output) {
  ::std::string ret;
  for (size_t at = 0; ; ) {
    const size_t line_end = output.find('\n', at);
    ret += "[  DEATH   ] ";
    if (line_end == ::std::string::npos) {
      ret += output.substr(at);
      break;
    }
    ret += output.substr(at, line_end + 1 - at);
    at = line_end + 1;
  }
  return ret;
}

#define GTEST_TEST_THROW_WHAT_(statement, expected_exception, re_pattern, fail) \
  GTEST_AMBIGUOUS_ELSE_BLOCKER_ \
  if (::testing::internal::ConstCharPtr gtest_msg = "") { \
    bool gtest_caught_expected = false; \
    try { \
      GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement); \
    } \
    catch (expected_exception const& exc) { \
      gtest_caught_expected = true; \
      char const * msg = exc.what();      \
      const ::testing::internal::RE& gtest_regex = (re_pattern); \
      if (msg == NULL || ! ::testing::internal::RE::PartialMatch(msg, gtest_regex ) ) { \
          std::ostringstream buffer; \
          buffer  << "    Result: throws but not with expected message.\n" \
                  << "  Expected: "  << gtest_regex.pattern() << "\n" \
                  << "Actual msg:\n" << FormatDeathTestOutput(std::string(msg)); \
          buffer << std::endl; \
          static std::string msg_buffer = buffer.str(); \
          gtest_msg.value = msg_buffer.c_str(); \
          goto GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__); \
      } \
    } \
    catch (...) { \
      gtest_msg.value = \
          "Expected: " #statement " throws an exception of type " \
          #expected_exception ".\n  Actual: it throws a different type."; \
      goto GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__); \
    } \
    if (!gtest_caught_expected) { \
      gtest_msg.value = \
          "Expected: " #statement " throws an exception of type " \
          #expected_exception ".\n  Actual: it throws nothing."; \
      goto GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__); \
    } \
  } else \
    GTEST_CONCAT_TOKEN_(gtest_label_testthrow_, __LINE__): \
      fail(gtest_msg.value)

#define EXPECT_THROW_WHAT(statement, expected_exception, pattern) \
  GTEST_TEST_THROW_WHAT_(statement, expected_exception, pattern,  GTEST_NONFATAL_FAILURE_)




#endif /* GTEST_THROW_WHAT_HH_ */
