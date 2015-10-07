/*
 * reader_to_storage_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#include <flow_gtest.hh>
#include <fstream>


#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "input/path_base.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"

using namespace std;

using namespace Input;

TEST(PathJSON, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());
    PathJSON path(in_str);

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
    EXPECT_EQ(1,path.find_ref_node()->get_int_value() );

    path.go_to_root();
    path.down(6);
    path.down("b");
    EXPECT_EQ("ctyri",path.find_ref_node()->get_string_value() );
}


TEST(PathJSON, errors) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";
	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());

    PathJSON path(in_str);
    path.down(8);
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcRefOfWrongType,"has wrong type, should by string." );

    path.go_to_root();
    path.down(9); // "REF":"/5/10"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "index out of size of Array" );

    path.go_to_root();
    path.down(10); // "REF":"/6/../.."
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "can not go up from root" );

    path.go_to_root();
    path.down(11); // "REF":"/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "there should be Record" );

    path.go_to_root();
    path.down(12); // "REF":"/6/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "key 'key' not found" );
}



TEST(PathYAML, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

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

    path.up();
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("b");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/b",os.str());
    }

}


TEST(PathYAML, values) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

	path.down(0); // bool value
	EXPECT_FALSE(path.get_bool_value());
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_STREQ("false", path.get_string_value().c_str());
	path.up();

	path.down(1); // int value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_EQ(1, path.get_int_value());
	EXPECT_FLOAT_EQ(1.0, path.get_double_value());
	EXPECT_STREQ("1", path.get_string_value().c_str());
	path.up();

	path.down(3); // double value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_FLOAT_EQ(3.3, path.get_double_value());
	EXPECT_STREQ("3.3", path.get_string_value().c_str());
	path.up();

	path.down(4); // string value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderToStorage::ExcInputError );
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
	path.up();

	path.down(9); // int64 value
	EXPECT_EQ(5000000000000, path.get_int_value());
	path.up();

	path.down(6); // record
	std::set<std::string> set;
	path.get_record_key_set(set);
	EXPECT_EQ(2, set.size());
	EXPECT_TRUE( set.find("a")!=set.end() );
	EXPECT_TRUE( set.find("b")!=set.end() );
	EXPECT_FALSE( set.find("c")!=set.end() );

	path.down("b"); // reference
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
}


class InputReaderToStorageTest : public testing::Test, public Input::ReaderToStorage {
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format = FileFormat::format_JSON) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	ReaderToStorage::read_stream(in, root_type, format);
    }
};

TEST_F(InputReaderToStorageTest, Integer) {
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

TEST_F(InputReaderToStorageTest, Double) {
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

TEST_F(InputReaderToStorageTest, Selection) {
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

TEST_F(InputReaderToStorageTest, String) {
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

TEST_F(InputReaderToStorageTest, Bool) {
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


const string input_yaml_array = R"YAML(
- 3.2
- 4
- 4.01
)YAML";


TEST_F(InputReaderToStorageTest, Array) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    Type::Array darr_type( Type::Double(3.1,4.1), 2,4);

    {  //JSON format
        stringstream ss("[ 3.2, 4, 4.01 ]");
        read_stream(ss, darr_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ(3.2, storage_->get_item(0)->get_double() );
        EXPECT_EQ(4, storage_->get_item(1)->get_double() );
    }

    {  //YAML format
        stringstream ss( input_yaml_array );
        read_stream(ss, darr_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ(3.2, storage_->get_item(0)->get_double() );
        EXPECT_EQ(4, storage_->get_item(1)->get_double() );
    }

    {  //JSON format
        stringstream ss("[ 3.2 ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ExcInputError, "Do not fit into size limits of the Array.");
    }

    {  //YAML format
        stringstream ss("- 3.2");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type, FileFormat::format_YAML);} , ExcInputError, "Do not fit into size limits of the Array.");
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


TEST_F(InputReaderToStorageTest, Record) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    static Type::Record rec_type = Type::Record( "SomeRec","desc.")
    	.declare_key("int_key", Type::Integer(0,5), Type::Default::obligatory(), "")
    	.declare_key("str_key", Type::String(), "")
		.close();

    { // JSON format
        stringstream ss("{ int_key=5 }");
        read_stream(ss, rec_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_EQ(true, storage_->get_item(1)->is_null() );
    }

    { // YAML format
        stringstream ss("int_key: 5");
        read_stream(ss, rec_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_EQ(true, storage_->get_item(1)->is_null() );
    }

    { // JSON format
        stringstream ss("{ str_key=\"ahoj\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ExcInputError, "Missing obligatory key 'int_key'.");
    }

    { // YAML format
        stringstream ss("str_key: ahoj}");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type, FileFormat::format_YAML);} , ExcInputError, "Missing obligatory key 'int_key'.");
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


const string input_yaml_abstract = R"YAML(
!EqDarcy
b_val: 10
a_val: prime
mesh: some.msh
)YAML";


TEST_F(InputReaderToStorageTest, AbstractRec) {
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

    EXPECT_FALSE(b_rec.is_finished());

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
        EXPECT_EQ("EqDarcy", storage_->get_item(0)->get_string());
        EXPECT_EQ("some.msh", storage_->get_item(1)->get_string() );
        EXPECT_EQ("prime", storage_->get_item(2)->get_string() );
        EXPECT_EQ(10, storage_->get_item(3)->get_int() );
    }

    {   //Try other correct type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=5.5, mesh=\"some.msh\" }");
        read_stream(ss, a_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ("EqTransp", storage_->get_item(0)->get_string());
        EXPECT_EQ("some.msh", storage_->get_item(1)->get_string() );
        EXPECT_EQ(5.5, storage_->get_item(2)->get_double() );
        EXPECT_EQ(4, storage_->get_item(3)->get_int() );
    }

    {   //Try YAML correct type
        stringstream ss(input_yaml_abstract);
        read_stream(ss, a_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ("EqDarcy", storage_->get_item(0)->get_string());
        EXPECT_EQ("some.msh", storage_->get_item(1)->get_string() );
        EXPECT_EQ("prime", storage_->get_item(2)->get_string() );
        EXPECT_EQ(10, storage_->get_item(3)->get_int() );
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

       { // YAML format
		   stringstream ss("20");
		   this->read_stream(ss, ar, FileFormat::format_YAML);

		   EXPECT_NE((void *)NULL, storage_);
		   storage_->get_array_size();
		   EXPECT_EQ(3, storage_->get_array_size());
		   EXPECT_EQ("BR", storage_->get_item(0)->get_string());
		   EXPECT_EQ(10, storage_->get_item(1)->get_int());
		   EXPECT_EQ(20, storage_->get_item(2)->get_int());
       }

       { // JSON format
		   stringstream ss("20");
		   this->read_stream(ss, ar);

		   EXPECT_NE((void *)NULL, storage_);
		   storage_->get_array_size();
		   EXPECT_EQ(3, storage_->get_array_size());
		   EXPECT_EQ("BR", storage_->get_item(0)->get_string());
		   EXPECT_EQ(10, storage_->get_item(1)->get_int());
		   EXPECT_EQ(20, storage_->get_item(2)->get_int());
       }
    }


}


const string input_json_multiple_inheritance = R"JSON(
{
  primary = { TYPE="Desc_B", b_val=1 },
  secondary = { TYPE="Desc_B", b_val=5 }
}
)JSON";


TEST_F(InputReaderToStorageTest, AbstractMultipleInheritance) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	Type::AbstractRecord a_rec1 = Type::AbstractRecord("Base1", "Base of equation records.").close();
	Type::AbstractRecord a_rec2 = Type::AbstractRecord("Base2", "Other base of equation records.").close();

	Type::Record rec_a = Type::Record("Desc_A", "First descendant")
			.derive_from(a_rec1)
			.declare_key("a_val", Type::String(), Type::Default::obligatory(), "")
			.close();
	Type::Record rec_b = Type::Record("Desc_B", "Second descendant")
			.derive_from(a_rec1)
			.derive_from(a_rec2)
			.declare_key("b_val", Type::Integer(), Type::Default::obligatory(), "")
			.close();
	Type::Record rec_c = Type::Record("Desc_C", "Third descendant")
			.derive_from(a_rec2)
			.declare_key("c_val", Type::Double(), Type::Default::obligatory(), "")
			.close();

	Type::Record root = Type::Record("problem", "Root record")
			.declare_key("primary", a_rec1, Type::Default::obligatory(), "")
			.declare_key("secondary", a_rec2, Type::Default::obligatory(), "")
			.close();

    {   // Try correct type
        stringstream ss(input_json_multiple_inheritance);
        read_stream(ss, root);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());

        EXPECT_EQ(2, storage_->get_item(0)->get_array_size());
        EXPECT_EQ("Desc_B", storage_->get_item(0)->get_item(0)->get_string());
        EXPECT_EQ(1, storage_->get_item(0)->get_item(1)->get_int());

        EXPECT_EQ(2, storage_->get_item(1)->get_array_size());
        EXPECT_EQ("Desc_B", storage_->get_item(1)->get_item(0)->get_string());
        EXPECT_EQ(5, storage_->get_item(1)->get_item(1)->get_int());
    }

}


/*TEST_F(InputReaderToStorageTest, AdHocAbstractRec) {
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

TEST(InputReaderToStorageTest_external, get_root_interface) {
    static Type::Record one_rec = Type::Record("One","")
    	.declare_key("one",Input::Type::Integer(),"")
		.close();
    one_rec.finish();

    ReaderToStorage json_reader("{ one=1 }", one_rec, FileFormat::format_JSON);
    Input::Record rec=json_reader.get_root_interface<Input::Record>();
    EXPECT_EQ(1, *(rec.find<int>("one")) );
    //json_reader.get_storage()->print(cout);

}

TEST_F(InputReaderToStorageTest, default_values) {
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


