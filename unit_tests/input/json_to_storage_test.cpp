/*
 * json_to_storage.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#include <flow_gtest.hh>
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

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, const Type::TypeBase &root_type) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	JSONToStorage::read_stream(in,root_type);
    }
};

TEST_F(InputJSONToStorageTest, Integer) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Integer int_type(1,10);
    Type::Integer any_int;

    {
        stringstream ss("5");
        read_stream(ss, int_type);

        EXPECT_EQ(5, storage_->get_int());
    }
    {
        stringstream ss("5000000000");
        EXPECT_THROW_WHAT( {read_stream(ss, any_int);} , ExcInputError, "Value out of bounds.");
    }
    {
        stringstream ss("0");
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ExcInputError, "Value out of bounds.");
    }
    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ExcInputError, "The value should be 'JSON int', but we found.* 'JSON object'");
    }
}

TEST_F(InputJSONToStorageTest, Double) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Double dbl_type(1.1,10.1);
    Type::Double any_double;

    {
        stringstream ss("5.5");
        read_stream(ss, dbl_type);

        EXPECT_EQ(5.5, storage_->get_double());
    }

    {
        stringstream ss("5000000000000");
        read_stream(ss, any_double);

        EXPECT_EQ(5e12, storage_->get_double());
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
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON object'");
    }
}

TEST_F(InputJSONToStorageTest, Selection) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Selection sel_type = Type::Selection("IntSelection")
    	.add_value(10,"ten","")
    	.add_value(1,"one","")
		.close();
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
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ExcInputError, "The value should be 'JSON string', but we found:.* 'JSON object'");
    }
}

TEST_F(InputJSONToStorageTest, String) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::String str_type;

    {
        stringstream ss("\"Important message\"");
        read_stream(ss, str_type);

        EXPECT_EQ("Important message", storage_->get_string());
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, str_type);} , ExcInputError, "The value should be 'JSON string', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, bool_type);} , ExcInputError, "The value should be 'JSON bool', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "The value should be 'JSON array', but we found:.* 'JSON object'");
    }

    {
        stringstream ss("[ 3.2, {} ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss1, darr_type);}, ExcInputError , "The value should be 'JSON real', but we found:.* 'JSON object'");
    }

    // test auto conversion failed
    {
        stringstream ss("3.2");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);}, ExcInputError , "Automatic conversion to array not allowed. The value should be 'JSON array', but we found:.* 'JSON real'");
    }
}


TEST_F(InputJSONToStorageTest, Record) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    static Type::Record rec_type = Type::Record( "SomeRec","desc.")
    	.declare_key("int_key", Type::Integer(0,5), Type::Default::obligatory(), "")
    	.declare_key("str_key", Type::String(), "")
		.close();

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
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "The value should be 'JSON object', but we found:.* 'JSON array'");
    }

    {
        stringstream ss("{ int_key=6 }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Value out of bounds.");
    }

    // test auto conversion
    {

        static Type::Record sub_rec = Type::Record( "SubRecord", "")
        	.declare_key("bool_key", Type::Bool(), Type::Default("false"), "")
        	.declare_key("int_key", Type::Integer(),  "")
        	.allow_auto_conversion("int_key")
			.close();
        sub_rec.finish();

        stringstream ss("123");
        read_stream(ss, sub_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_FALSE( storage_->get_item(0)->get_bool() );
        EXPECT_EQ(123, storage_->get_item(1)->get_int() );

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ExcAutomaticConversionError , "The value should be 'JSON int', but we found:.* 'JSON real'");
    }

    {
        static Type::Record sub_rec = Type::Record( "SubRecord", "")
        	.close();
        sub_rec.finish();

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ExcInputError , "The value should be 'JSON object', but we found:.* 'JSON real'");
    }

    // Test automatic conversion from record
/*
    {
        static Type::Record lower( "Lower", "");
        lower.declare_key("int", Type::Integer(), "");
        lower.finish();

        static Type::Record upper( "Upper", "");
        upper.has_obligatory_type_key();
        upper.declare_key("rec", lower, "");
        upper.allow_auto_conversion("rec");
        upper.finish();

        stringstream ss("{TYPE=\"Upper\", rec={int=37} }");
        read_stream(ss, upper);
        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ( 0, storage_->get_item(0)->get_int() );
        EXPECT_EQ( 37, storage_->get_item(1)->get_item(0)->get_int() );

        stringstream ss1("{int=37}");
        read_stream(ss1, upper);
        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ( 0, storage_->get_item(0)->get_int() );
        EXPECT_EQ( 37, storage_->get_item(1)->get_item(0)->get_int() );

    }
*/
/*
    {
        static Type::AbstractRecord abstr("Abstract", "");
        abstr.finish();

        static Type::Record lower( "Lower", "");
        lower.derive_from(abstr);
        lower.declare_key("int", Type::Integer(), "");
        lower.finish();

        static Type::Record upper( "Upper", "");
        upper.decalare_key("rec", upper, "");
        upper.allow_auto_conversion("rec");
        upper.has_obligatory_type_key();
        upper.finish();

        stringstream ss("{TYPE=\"Upper\", rec={int=37} }");
        read_stream(ss, upper);
        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ( 0, storage_->get_item(0)->get_int() );
        EXPECT_EQ( 37, storage_->get_item(1)->get_item(0) );

        stringstream ss("{int=37}");
        read_stream(ss, upper);
        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ( 0, storage_->get_item(0)->get_int() );
        EXPECT_EQ( 37, storage_->get_item(1)->get_item(0) );

    }
*/

}

TEST_F(InputJSONToStorageTest, AbstractRec) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

	static Type::Record copy_rec = Type::Record("Copy","")
       	.declare_key("mesh", Type::String(), Type::Default("input.msh"), "Comp. mesh.")
       	.declare_key("a_val", Type::String(), Type::Default::obligatory(), "")
		.close();

    static Type::AbstractRecord a_rec = Type::AbstractRecord("EqBase","Base of equation records.")
    	.close();

    static Type::Record b_rec = Type::Record("EqDarcy","")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
    	.declare_key("b_val", Type::Integer(), "")
		.close();

    EXPECT_EQ(false, b_rec.is_finished());

    static Type::Record c_rec = Type::Record("EqTransp","")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
		.declare_key("c_val", Type::Integer(), "")
		.declare_key("a_val", Type::Double(),"")
		.close();

    a_rec.add_child(b_rec);
    a_rec.add_child(c_rec);
    a_rec.finish();
    b_rec.finish();
    c_rec.finish();

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
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON string'");

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
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ExcInputError, "The value should be 'JSON object', but we found:.* 'JSON array'");

    }

    { // auto conversion
       Type::AbstractRecord ar = Type::AbstractRecord("AR","")
         .allow_auto_conversion("BR")
         .close();
       Type::Record br = Type::Record("BR","")
         .derive_from(ar)
         .declare_key("x",Input::Type::Integer(),Input::Type::Default("10"),"")
         .declare_key("y",Input::Type::Integer(),"")
         .allow_auto_conversion("y")
         .close();
       ar.add_child(br);
       br.finish();

       stringstream ss("20");
       this->read_stream(ss, ar);

       EXPECT_NE((void *)NULL, storage_);
       storage_->get_array_size();
       EXPECT_EQ(3, storage_->get_array_size());
       EXPECT_EQ(0, storage_->get_item(0)->get_int());
       EXPECT_EQ(10, storage_->get_item(1)->get_int());
       EXPECT_EQ(20, storage_->get_item(2)->get_int());
    }


}

/*TEST_F(InputJSONToStorageTest, AdHocAbstractRec) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Selection sel_type("TYPE_selection");
    sel_type.add_value(0, "EqDarcy");
    sel_type.add_value(1, "EqTransp");
    sel_type.finish();

    static Type::Record b_rec("EqDarcy","");
    b_rec.declare_key("TYPE", sel_type, Type::Default("EqDarcy"), "Type of problem");
    b_rec.declare_key("a_val", Type::String(), Type::Default("Description"), "");
    b_rec.declare_key("b_val", Type::Integer(), "");
    b_rec.declare_key("mesh", Type::String(), Type::Default::obligatory(), "Mesh.");
    b_rec.finish();

    static Type::Record c_rec("EqTransp","");
    c_rec.declare_key("TYPE", sel_type, Type::Default("EqTransp"), "Type of problem");
    c_rec.declare_key("a_val", Type::Double(),"");
    c_rec.declare_key("c_val", Type::Integer(), "");
    c_rec.finish();

    static Type::AbstractRecord a_rec("EqBase","Base of equation records.");
    static Type::AdHocAbstractRecord ah_rec(a_rec);
    ah_rec.add_child(b_rec);
    ah_rec.add_child(c_rec);
    ah_rec.finish();

    {   // Try one correct type
        stringstream ss("{ TYPE=\"EqDarcy\", b_val=4, mesh=\"some.msh\" }");
        read_stream(ss, ah_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(0, storage_->get_item(0)->get_int());
        EXPECT_EQ("Description", storage_->get_item(1)->get_string() );
        EXPECT_EQ(4, storage_->get_item(2)->get_int() );
        EXPECT_EQ("some.msh", storage_->get_item(3)->get_string() );
    }

    {   //Try other correct type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=5.5 }");
        read_stream(ss, ah_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ(1, storage_->get_item(0)->get_int());
        EXPECT_EQ(5.5, storage_->get_item(1)->get_double() );
        EXPECT_EQ(4, storage_->get_item(2)->get_int() );
    }

    {   // Missing TYPE
        stringstream ss("{ b_val=4, a_val=\"Some text\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, ah_rec);}, ExcInputError, "Missing key 'TYPE' in AbstractRecord.");
    }

    {   // Wrong derived value type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=\"prime\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, ah_rec);}, ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON string'");
    }
} */

TEST(InputJSONToStorageTest_external, get_root_interface) {
    static Type::Record one_rec = Type::Record("One","")
    	.declare_key("one",Input::Type::Integer(),"")
		.close();
    one_rec.finish();

    stringstream ss("{ one=1 }");
    JSONToStorage json_reader(ss, one_rec);
    Input::Record rec=json_reader.get_root_interface<Input::Record>();
    EXPECT_EQ(1, *(rec.find<int>("one")) );
    //json_reader.get_storage()->print(cout);

}

TEST_F(InputJSONToStorageTest, default_values) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Selection sel_type = Type::Selection("tmp selection")
    	.add_value(1,"one")
		.add_value(2,"two")
		.close();
    sel_type.finish();

    static Type::Record sub_rec = Type::Record( "SubRecord", "")
    	.declare_key("bool_key", Type::Bool(), Type::Default("false"), "")
    	.declare_key("int_key", Type::Integer(),  "")
    	.allow_auto_conversion("int_key")
    	.close();

    static Type::Record rec_type = Type::Record( "SomeRec","desc.")
    	.declare_key("int_key", Type::Integer(0,5), Type::Default("4"), "")
    	.declare_key("bool_key", Type::Bool(), Type::Default("true"),"")
    	.declare_key("sel_key", sel_type, Type::Default("two"),"")
    	.declare_key("double_key", Type::Double(), Type::Default("1.23"),"")
    	.declare_key("str_key", Type::String(), Type::Default("ahoj"),"")
    	.declare_key("array_key", Type::Array( Type::Integer() ), Type::Default("123"), "")
	    .declare_key("rec_key", sub_rec, Type::Default("321"), "")
		.close();

    rec_type.finish();
    sub_rec.finish();

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


