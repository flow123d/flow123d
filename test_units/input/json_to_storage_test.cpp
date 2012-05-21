/*
 * json_to_storage.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test method get_root_interface
 * TODO: move EXPECT_THROW_WHAT_ into flow_test.hh
 * TODO: test catching of errors in JSON file format.
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
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ExcInputError, "Value out of bounds.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ExcInputError, "Wrong type, has to be Int.");
    }
}

TEST_F(InputJSONToStorageTest, Double) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Double dbl_type(1.1,10.1);

    {
        stringstream ss("5.5");
        read_stream(ss, dbl_type);

        EXPECT_EQ(5.5, storage_->get_double());
    }

    {
        stringstream ss("5");
        read_stream(ss, dbl_type);

        EXPECT_EQ(5, storage_->get_double());
    }

    {
        stringstream ss("0");
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ExcInputError, "Value out of bounds.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ExcInputError, "Wrong type, has to be Double.");
    }
}

TEST_F(InputJSONToStorageTest, Selection) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Selection<int> sel_type("IntSelection");
    sel_type.add_value(10,"ten","");
    sel_type.add_value(1,"one","");
    sel_type.finish();

    {
        stringstream ss("\"ten\"");
        read_stream(ss, sel_type);

        EXPECT_EQ(10, storage_->get_int());
    }

    {
        stringstream ss("\"red\"");
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ExcInputError, "Wrong value of the Selection.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ExcInputError, "Wrong type, value should be String .key of Selection.");
    }
}

TEST_F(InputJSONToStorageTest, String) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::String str_type;

    {
        F_ENTRY;
        stringstream ss("\"Important message\"");
        read_stream(ss, str_type);

        EXPECT_EQ("Important message", storage_->get_string());
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, str_type);} , ExcInputError, "Wrong type, has to be String.");
    }
}

TEST_F(InputJSONToStorageTest, Bool) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Bool bool_type;

    {
        stringstream ss("true");
        read_stream(ss, bool_type);

        EXPECT_EQ(true, storage_->get_bool());
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, bool_type);} , ExcInputError, "Wrong type, has to be Bool.");
    }
}

TEST_F(InputJSONToStorageTest, Array) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Array darr_type( Type::Double(3.1,4.1), 2,4);

    {
        stringstream ss("[ 3.2, 4, 4.01 ]");
        read_stream(ss, darr_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ(3.2, storage_->get_item(0)->get_double() );
        EXPECT_EQ(4, storage_->get_item(1)->get_double() );
    }

    {
        stringstream ss("[ 3.2 ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Do not fit into size limits of the Array.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Wrong type, has to be Array.");
    }

    {
        stringstream ss("[ 3.2, {} ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Wrong type, has to be Double.");
    }

    {
        stringstream ss("[ 3.0, 3.9 ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Value out of bounds.");
    }

}


TEST_F(InputJSONToStorageTest, Record) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    static Type::Record rec_type( "SomeRec","desc.");
    rec_type.declare_key("int_key", Type::Integer(0,5), Type::DefaultValue(Type::DefaultValue::obligatory), "");
    rec_type.declare_key("str_key", Type::String(), "");
    rec_type.finish();

    {
        stringstream ss("{ int_key=5 }");
        read_stream(ss, rec_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_EQ(true, storage_->get_item(1)->is_null() );
    }

    {
        stringstream ss("{ str_key=\"ahoj\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Missing obligatory key 'int_key'.");
    }

    {
        stringstream ss("[]");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Wrong type, has to be Record.");
    }

    {
        stringstream ss("{ int_key=6 }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Value out of bounds.");
    }
}

TEST_F(InputJSONToStorageTest, AbstratRec) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::AbstractRecord a_rec("EqBase","Base of equation records.");
    a_rec.declare_key("mesh", Type::String(), Type::DefaultValue(Type::DefaultValue::obligatory), "Comp. mesh.");
    a_rec.declare_key("a_val", Type::String(), Type::DefaultValue(Type::DefaultValue::obligatory), "");
    a_rec.finish();

    static Type::Record b_rec("EqDarcy","");
    b_rec.derive_from(a_rec);
    b_rec.declare_key("b_val", Type::Integer(), "");
    b_rec.finish();

    EXPECT_EQ(true, b_rec.is_finished());

    static Type::Record c_rec("EqTransp","");
    c_rec.derive_from(a_rec);
    c_rec.declare_key("c_val", Type::Integer(), "");
    c_rec.declare_key("a_val", Type::Double(),"");
    c_rec.finish();

    a_rec.no_more_descendants();
    EXPECT_EQ(true, b_rec.is_finished());
    EXPECT_EQ(b_rec, a_rec.get_descendant("EqDarcy"));
    EXPECT_EQ(true, a_rec.get_descendant("EqDarcy").is_finished());
    EXPECT_EQ(Type::Double() , *( c_rec.key_iterator("a_val")->type_));
    //cout << a_rec;


    {   // Try one correct type
        stringstream ss("{ TYPE=\"EqDarcy\", b_val=10, a_val=\"prime\", mesh=\"some.msh\" }");
        read_stream(ss, a_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(0, storage_->get_item(0)->get_int());
        EXPECT_EQ("some.msh", storage_->get_item(1)->get_string() );
        EXPECT_EQ("prime", storage_->get_item(2)->get_string() );
        EXPECT_EQ(10, storage_->get_item(3)->get_int() );
    }

    {   //Try other correct type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=5.5, mesh=\"some.msh\" }");
        read_stream(ss, a_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(1, storage_->get_item(0)->get_int());
        EXPECT_EQ("some.msh", storage_->get_item(1)->get_string() );
        EXPECT_EQ(5.5, storage_->get_item(2)->get_double() );
        EXPECT_EQ(4, storage_->get_item(3)->get_int() );
    }

    {   // Wrong derived value type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Wrong type, has to be Double.");

    }

    {   // Missing TYPE
        stringstream ss("{ c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Missing key 'TYPE' in AbstractRecord.");

    }

    {   // Wrong derived value type
        stringstream ss("{ TYPE=\"EqTrans\", c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Wrong TYPE='EqTrans' of AbstractRecord.");

    }

    {   // Wrong derived value type
        stringstream ss("[]");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Wrong type, has to be Record.");

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



