/*
 * json_to_storage.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#include <gtest/gtest.h>
#include <gtest_throw_what.hh>
#include <fstream>

#include "input/json_to_storage.hh"

using namespace std;

using namespace Input;

TEST(JSONPath, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/json_to_storage_test.con").c_str());
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


TEST(JSONPath, errors) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";
	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/json_to_storage_test.con").c_str());
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
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ExcInputError, "The value should be 'JSON int', but we found type: 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ExcInputError, "The value should be 'JSON real', but we found type: 'JSON object'");
    }
}

TEST_F(InputJSONToStorageTest, Selection) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Selection sel_type("IntSelection");
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
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ExcInputError, "Wrong value 'red' of the Selection.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ExcInputError, "The value should be 'JSON string', but we found type: 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, str_type);} , ExcInputError, "The value should be 'JSON string', but we found type: 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, bool_type);} , ExcInputError, "The value should be 'JSON bool', but we found type: 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "The value should be 'JSON array', but we found type: 'JSON object'");
    }

    {
        stringstream ss("[ 3.2, {} ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "The value should be 'JSON real', but we found type: 'JSON object'");
    }

    {
        stringstream ss("[ 3.0, 3.9 ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Value out of bounds.");
    }

    // test automatic conversion
    {
        Type::Array darr_type( Type::Double(3.1,4.1));
        stringstream ss("3.2");
        read_stream(ss, darr_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(1, storage_->get_array_size());
        EXPECT_EQ(3.2, storage_->get_item(0)->get_double() );

        stringstream ss1("{ key=3.2}");
        EXPECT_THROW_WHAT( {read_stream(ss1, darr_type);}, ExcInputError , "The value should be 'JSON real', but we found type: 'JSON object'");
    }

    // test auto conversion failed
    {
        stringstream ss("3.2");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);}, ExcInputError , "Automatic conversion to array not allowed. The value should be 'JSON array', but we found type: 'JSON real'");
    }
}


TEST_F(InputJSONToStorageTest, Record) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    static Type::Record rec_type( "SomeRec","desc.");
    rec_type.declare_key("int_key", Type::Integer(0,5), Type::Default::obligatory(), "");
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
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "The value should be 'JSON object', but we found type: 'JSON array'");
    }

    {
        stringstream ss("{ int_key=6 }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Value out of bounds.");
    }

    // test auto conversion
    {

        static Type::Record sub_rec( "SubRecord", "");
        sub_rec.declare_key("bool_key", Type::Bool(), Type::Default("false"), "");
        sub_rec.declare_key("int_key", Type::Integer(),  "");
        sub_rec.allow_auto_conversion("int_key");
        sub_rec.finish();

        stringstream ss("123");
        read_stream(ss, sub_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_FALSE( storage_->get_item(0)->get_bool() );
        EXPECT_EQ(123, storage_->get_item(1)->get_int() );

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ExcInputError , "The value should be 'JSON int', but we found type: 'JSON real'");
    }

    {
        static Type::Record sub_rec( "SubRecord", "");
        sub_rec.finish();

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ExcInputError , "The value should be 'JSON object', but we found type: 'JSON real'");
    }


}

TEST_F(InputJSONToStorageTest, AbstractRec) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::AbstractRecord a_rec("EqBase","Base of equation records.");
    a_rec.declare_key("mesh", Type::String(), Type::Default::obligatory(), "Comp. mesh.");
    a_rec.declare_key("a_val", Type::String(), Type::Default::obligatory(), "");
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
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "The value should be 'JSON real', but we found type: 'JSON string'");

    }

    {   // Missing TYPE
        stringstream ss("{ c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Missing key 'TYPE' in AbstractRecord.");

    }

    {   // Wrong value of TYPE
        stringstream ss("{ TYPE=\"EqTrans\", c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "Wrong value 'EqTrans' of the Selection.");

    }

    {   // Wrong derived value type
        stringstream ss("[]");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "The value should be 'JSON object', but we found type: 'JSON array'");

    }


}

TEST(InputJSONToStorageTest_external, get_root_interface) {
    static Input::Type::Record one_rec("One","");
    one_rec.declare_key("one",Input::Type::Integer(),"");
    one_rec.finish();

    stringstream ss("{ one=1 }");
    JSONToStorage json_reader;
    json_reader.read_stream(ss, one_rec);
    Input::Record rec=json_reader.get_root_interface<Input::Record>();
    EXPECT_EQ(1, *(rec.find<int>("one")) );
    //json_reader.get_storage()->print(cout);

}

TEST_F(InputJSONToStorageTest, default_values) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Selection sel_type("tmp selection");
    sel_type.add_value(1,"one");
    sel_type.add_value(2,"two");
    sel_type.finish();

    static Type::Record rec_type( "SomeRec","desc.");
    rec_type.declare_key("int_key", Type::Integer(0,5), Type::Default("4"), "");
    rec_type.declare_key("bool_key", Type::Bool(), Type::Default("true"),"");
    rec_type.declare_key("sel_key", sel_type, Type::Default("two"),"");
    rec_type.declare_key("double_key", Type::Double(), Type::Default("1.23"),"");
    rec_type.declare_key("str_key", Type::String(), Type::Default("ahoj"),"");
    rec_type.declare_key("array_key", Type::Array( Type::Integer() ), Type::Default("123"), "");

    static Type::Record sub_rec( "SubRecord", "");
    sub_rec.declare_key("bool_key", Type::Bool(), Type::Default("false"), "");
    sub_rec.declare_key("int_key", Type::Integer(),  "");
    sub_rec.allow_auto_conversion("int_key");
    sub_rec.finish();

    rec_type.declare_key("rec_key", sub_rec, Type::Default("321"), "");
    rec_type.finish();

    {
        stringstream ss("{ }");
        read_stream(ss, rec_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(7, storage_->get_array_size());
        EXPECT_EQ(4, storage_->get_item(0)->get_int() );
        EXPECT_TRUE( storage_->get_item(1)->get_bool() );
        EXPECT_EQ(2, storage_->get_item(2)->get_int() );
        EXPECT_EQ(1.23, storage_->get_item(3)->get_double() );
        EXPECT_EQ("ahoj", storage_->get_item(4)->get_string() );
        EXPECT_EQ(123 , storage_->get_item(5)->get_item(0)->get_int() );
        EXPECT_FALSE( storage_->get_item(6)->get_item(0)->get_bool() );
        EXPECT_EQ(321 , storage_->get_item(6)->get_item(1)->get_int() );
    }
}


