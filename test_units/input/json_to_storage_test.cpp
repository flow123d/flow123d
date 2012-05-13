/*
 * json_to_storage.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#include <gtest/gtest.h>

#include "json_to_storage.hh"

using namespace std;

// Test data strings
const string storage_input = R"JSON(
[ 
  false,
  1,
  2,
  3.3,
  "ctyri",
  [1,2,3,4,5],
  { "a": {"REF":"../../7/b/c"},
    "b": {"REF":"/4"}
  },
  { "b": { "c":1 } },
  {"REF":[0]},          
  {}
]
)JSON";

using namespace Input;

TEST(JSONPath, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    stringstream in_str(storage_input);
    JSONPath::Node node;
    json_spirit::read_or_throw( in_str, node);
    JSONPath path(node);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down(6);
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("a");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }

    string ref;
    path.get_ref_from_head(ref);
    EXPECT_EQ("../../7/b/c",ref);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }

    EXPECT_EQ(1,path.find_ref_node(ref).head()->get_int() );

    path=path.find_ref_node("/6/b");
    path.get_ref_from_head(ref);
    EXPECT_EQ("/4", ref);
    EXPECT_EQ("ctyri",path.find_ref_node(ref).head()->get_str() );
}

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

TEST(JSONPath, errors) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";
    stringstream in_str(storage_input);
    JSONPath::Node node;
    json_spirit::read_or_throw( in_str, node);

    string ref;
    JSONPath path(node);
    path.down(8);

    EXPECT_THROW_WHAT( { path.get_ref_from_head(ref);} ,JSONPath::ExcRefOfWrongType,"has wrong type, should by string." );
    EXPECT_THROW_WHAT( { path.find_ref_node("/6/4");}, JSONPath::ExcReferenceNotFound, "there should be Array" );
    EXPECT_THROW_WHAT( { path.find_ref_node("/5/10");}, JSONPath::ExcReferenceNotFound, "index out of size of Array" );
    EXPECT_THROW_WHAT( { path.find_ref_node("/6/../..");}, JSONPath::ExcReferenceNotFound, "can not go up from root" );
    EXPECT_THROW_WHAT( { path.find_ref_node("/key");}, JSONPath::ExcReferenceNotFound, "there should be Record" );
    EXPECT_THROW_WHAT( { path.find_ref_node("/6/key");}, JSONPath::ExcReferenceNotFound, "key 'key' not found" );
}


class InputJSONToStorageTest : public testing::Test, public Input::JSONToStorage {

protected:

    virtual void SetUp() {
    }

    virtual void TearDown() {
    };

    ::Input::Type::Record *main;
    ::Input::Interface::Storage * storage;
};

TEST_F(InputJSONToStorageTest, Integer) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Integer int_type(1,10);

    {
        stringstream ss("5");
        read_stream(ss, int_type);
        EXPECT_EQ(5, storage_->get_int());
    }

    {
        stringstream ss("0");
        EXPECT_DEATH( {read_stream(ss, int_type);} , "out of bounds.");
    }

    {
        stringstream ss("{}");
        EXPECT_DEATH( {read_stream(ss, int_type);} , "has to be of type Int.");
    }
}
/*
    stringstream in_str(storage_input);
    Input::Type::Record rec("SomeRecord", "");

    Input::Storage storage(in_str, rec);

    EXPECT_FALSE( storage.get_item(0).get_bool());
    EXPECT_EQ(1, storage.get_item(1).get_int());
    EXPECT_EQ(2, storage.get_item(2).get_int());
    EXPECT_EQ(3.3, storage.get_item(3).get_double());
    EXPECT_EQ("ctyri", storage.get_item(4).get_string());
    EXPECT_EQ(5, storage.get_item(5).get_array_size());
    EXPECT_EQ(4, storage.get_item(5).get_item(3).get_int());
*/



