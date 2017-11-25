/*
 * reader_to_storage_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <fstream>


#include "input/reader_to_storage.hh"
#include "input/reader_internal_base.hh"
#include "input/accessors.hh"

using namespace std;

using namespace Input;


class InputReaderToStorageTest : public testing::Test, public Input::ReaderToStorage {
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

    // overload parent class method in order to reset pointers
    void read_stream(istream &in, Type::TypeBase &root_type, FileFormat format = FileFormat::format_JSON) {
    	this->storage_ = nullptr;
    	this->root_type_ = nullptr;
    	root_type.finish();
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
        EXPECT_THROW_WHAT( {read_stream(ss, any_int);} , ReaderInternalBase::ExcInputError, "Value out of bounds.");
    }
    {
        stringstream ss("0");
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ReaderInternalBase::ExcInputError, "Value out of bounds.");
    }
    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, int_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON int', but we found.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ReaderInternalBase::ExcInputError, "Value out of bounds.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, dbl_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ReaderInternalBase::ExcInputError, "Wrong value 'red' of the Selection.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, sel_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON string', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, str_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON string', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, bool_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON bool', but we found:.* 'JSON object'");
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
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ReaderInternalBase::ExcInputError, "Do not fit the size 1 of the Array.");
    }

    {  //YAML format
        stringstream ss("- 3.2");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type, FileFormat::format_YAML);} , ReaderInternalBase::ExcInputError, "Do not fit the size 1 of the Array.");
    }

    {
        stringstream ss("{}");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON array', but we found:.* 'JSON object'");
    }

    {
        stringstream ss("[ 3.2, {} ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON object'");
    }

    {
        stringstream ss("[ 3.0, 3.9 ]");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);} , ReaderInternalBase::ExcInputError, "Value out of bounds.");
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
        EXPECT_THROW_WHAT( {read_stream(ss1, darr_type);}, ReaderInternalBase::ExcInputError , "The value should be 'JSON real', but we found:.* 'JSON object'");
    }

    // test auto conversion failed
    {
        stringstream ss("3.2");
        EXPECT_THROW_WHAT( {read_stream(ss, darr_type);}, ReaderInternalBase::ExcInputError , "conversion to the single element array not allowed. Require type: 'JSON array'.* 'JSON real'");
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
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ("SomeRec", storage_->get_item(0)->get_string() );
        EXPECT_EQ(5, storage_->get_item(1)->get_int() );
        EXPECT_EQ(true, storage_->get_item(2)->is_null() );
    }

    { // YAML format
        stringstream ss("int_key: 5");
        read_stream(ss, rec_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ("SomeRec", storage_->get_item(0)->get_string() );
        EXPECT_EQ(5, storage_->get_item(1)->get_int() );
        EXPECT_EQ(true, storage_->get_item(2)->is_null() );
    }

    { // JSON format
        stringstream ss("{ str_key=\"ahoj\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ReaderInternalBase::ExcInputError, "Missing obligatory key 'int_key'.");
    }

    { // YAML format
        stringstream ss("str_key: ahoj}");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type, FileFormat::format_YAML);} , ReaderInternalBase::ExcInputError, "Missing obligatory key 'int_key'.");
    }

    {
        stringstream ss("[]");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ReaderInternalBase::ExcInputError, "The value should be 'JSON object', but we found:.* 'JSON array'");
    }

    {
        stringstream ss("{ int_key=6 }");
        EXPECT_THROW_WHAT( {read_stream(ss, rec_type);} , ReaderInternalBase::ExcInputError, "Value out of bounds.");
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
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_FALSE( storage_->get_item(1)->get_bool() );
        EXPECT_EQ(123, storage_->get_item(2)->get_int() );

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ReaderInternalBase::ExcAutomaticConversionError , "The value should be 'JSON int', but we found:.* 'JSON real'");
    }

    {
        static Type::Record sub_rec = Type::Record( "SubRecord", "")
        	.close();
        sub_rec.finish();

        stringstream ss1("1.23");
        EXPECT_THROW_WHAT( {read_stream(ss1, sub_rec);}, ReaderInternalBase::ExcInputError , "The value should be 'JSON object', but we found:.* 'JSON real'");
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
        static Type::Abstract abstr("Abstract", "");
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
       	.declare_key("mesh", Type::String(), Type::Default("\"input.msh\""), "Comp. mesh.")
       	.declare_key("a_val", Type::String(), Type::Default::obligatory(), "")
		.close();

    static Type::Abstract a_rec = Type::Abstract("EqBase","Base of equation records.")
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
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ReaderInternalBase::ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON string'");

    }

    {   // Missing TYPE
        stringstream ss("{ c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ReaderInternalBase::ExcInputError, "Can not determine type of the Abstract.");

    }

    {   // Wrong value of TYPE
        stringstream ss("{ TYPE=\"EqTrans\", c_val=4, a_val=\"prime\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ReaderInternalBase::ExcInputError, "Wrong value 'EqTrans' of the Selection.");

    }

    {   // Wrong derived value type
        stringstream ss("[]");
        EXPECT_THROW_WHAT( {read_stream(ss, a_rec);}, ReaderInternalBase::ExcInputError, "Can not determine type of the Abstract.");

    }

    { // auto conversion
       Type::Abstract ar = Type::Abstract("AR","")
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

	Type::Abstract a_rec1 = Type::Abstract("Base1", "Base of equation records.").close();
	Type::Abstract a_rec2 = Type::Abstract("Base2", "Other base of equation records.").close();

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
        EXPECT_EQ(3, storage_->get_array_size());

        EXPECT_EQ("problem", storage_->get_item(0)->get_string());

        EXPECT_EQ(2, storage_->get_item(1)->get_array_size());
        EXPECT_EQ("Desc_B", storage_->get_item(1)->get_item(0)->get_string());
        EXPECT_EQ(1, storage_->get_item(1)->get_item(1)->get_int());

        EXPECT_EQ(2, storage_->get_item(2)->get_array_size());
        EXPECT_EQ("Desc_B", storage_->get_item(2)->get_item(0)->get_string());
        EXPECT_EQ(5, storage_->get_item(2)->get_item(1)->get_int());
    }

}


TEST_F(InputReaderToStorageTest, AdHocAbstract) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Abstract a_rec = Type::Abstract("EqBase","Base Abstract of equation records.")
    	.close();

    static Type::Record b_rec = Type::Record("EqDarcy","")
    	.declare_key("a_val", Type::String(), Type::Default("\"Description\""), "")
    	.declare_key("b_val", Type::Integer(), "")
    	.declare_key("mesh", Type::String(), Type::Default::obligatory(), "Mesh.")
		.close();

    static Type::Record c_rec = Type::Record("EqTransp","")
    	.declare_key("a_val", Type::Double(),"")
    	.declare_key("c_val", Type::Integer(), "")
    	.close();

    static Type::AdHocAbstract ah_rec(a_rec);
    ah_rec.add_child(b_rec);
    ah_rec.add_child(c_rec);
    ah_rec.finish();

    {   // Try one correct type
        stringstream ss("{ TYPE=\"EqDarcy\", b_val=4, mesh=\"some.msh\" }");
        read_stream(ss, ah_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ("EqDarcy", storage_->get_item(0)->get_string());
        EXPECT_EQ("Description", storage_->get_item(1)->get_string() );
        EXPECT_EQ(4, storage_->get_item(2)->get_int() );
        EXPECT_EQ("some.msh", storage_->get_item(3)->get_string() );
    }

    {   //Try other correct type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=5.5 }");
        read_stream(ss, ah_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3, storage_->get_array_size());
        EXPECT_EQ("EqTransp", storage_->get_item(0)->get_string());
        EXPECT_EQ(5.5, storage_->get_item(1)->get_double() );
        EXPECT_EQ(4, storage_->get_item(2)->get_int() );
    }

    {   // Missing TYPE
        stringstream ss("{ b_val=4, a_val=\"Some text\", mesh=\"some.msh\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, ah_rec);}, ReaderInternalBase::ExcInputError, "Can not determine type of the Abstract.");
    }

    {   // Wrong derived value type
        stringstream ss("{ TYPE=\"EqTransp\", c_val=4, a_val=\"prime\" }");
        EXPECT_THROW_WHAT( {read_stream(ss, ah_rec);}, ReaderInternalBase::ExcInputError, "The value should be 'JSON real', but we found:.* 'JSON string'");
    }
}


const string input_yaml_tuple = R"YAML(
- 5
- 2.0
- some string
)YAML";

const string input_yaml_tuple_as_rec = R"YAML(
int_key: 5
str_key: some string
dbl_key: 2.0
)YAML";

const string input_yaml_tuple_with_null = R"YAML(
- 5
- null
- some string
)YAML";

TEST_F(InputReaderToStorageTest, Tuple) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    static Type::Tuple tuple_type = Type::Tuple( "SomeTuple","desc.")
    	.declare_key("int_key", Type::Integer(0,10), Type::Default::obligatory(), "")
    	.declare_key("dbl_key", Type::Double(0.0), Type::Default("1.0"), "")
    	.declare_key("str_key", Type::String(), Type::Default::optional(), "")
    	.declare_key("bool_key", Type::Bool(), Type::Default::optional(), "")
		.close();

    { // YAML format, Tuple defined as array
        stringstream ss( input_yaml_tuple );
        read_stream(ss, tuple_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_FLOAT_EQ(2.0, storage_->get_item(1)->get_double() );
        EXPECT_STREQ("some string", storage_->get_item(2)->get_string().c_str() );
    }

    { // YAML format, Tuple defined as record
        stringstream ss( input_yaml_tuple_as_rec );
        read_stream(ss, tuple_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_FLOAT_EQ(2.0, storage_->get_item(1)->get_double() );
        EXPECT_STREQ("some string", storage_->get_item(2)->get_string().c_str() );
    }

    { // YAML format, Tuple defined as single value (used auto-conversion)
        stringstream ss("5");
        read_stream(ss, tuple_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_FLOAT_EQ(1.0, storage_->get_item(1)->get_double() );
    }

    { // YAML format, Tuple contains null value
        stringstream ss( input_yaml_tuple_with_null );
        read_stream(ss, tuple_type, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(4, storage_->get_array_size());
        EXPECT_EQ(5, storage_->get_item(0)->get_int() );
        EXPECT_FLOAT_EQ(1.0, storage_->get_item(1)->get_double() );
        EXPECT_STREQ("some string", storage_->get_item(2)->get_string().c_str() );
    }

    static Type::Tuple int_tuple_type = Type::Tuple( "NumericTuple","desc.")
    	.declare_key("int1_key", Type::Integer(0,10), Type::Default::obligatory(), "")
    	.declare_key("int2_key", Type::Integer(0,20), Type::Default::obligatory(), "")
    	.declare_key("int3_key", Type::Integer(), Type::Default::optional(), "")
		.close();

    { // YAML format, Tuple defined as single value throws exception
        stringstream ss("- 5");
        EXPECT_THROW_WHAT( { read_stream(ss, int_tuple_type, FileFormat::format_YAML); },
        		ReaderInternalBase::ExcInputError, "Tuple with 2 obligatory keys" );
    }

}



const string input_yaml_empty_rec_without_tag = R"YAML(
format: 
  #int_key: 10
)YAML";

const string input_yaml_empty_rec_with_tag = R"YAML(
format: !RecordA
  #int_key: 10
)YAML";

const string input_yaml_rec_with_empty_scalar = R"YAML(
format: 
  int_key:
  str_key:
)YAML";



TEST_F(InputReaderToStorageTest, EmptyRecord) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Record inrec = Type::Record( "InRecord", "desc.")
        .declare_key("int_key", Type::Integer(), Type::Default::optional(), "")
        .declare_key("str_key", Type::String(), Type::Default::optional(), "")
        .close();

    static Type::Record root_rec = Type::Record( "RootRecord", "desc.")
        .declare_key("format", inrec, Type::Default::obligatory(), "")
	    .close();

    static Type::Abstract a_rec = Type::Abstract("BaseOfRec", "Base of records.")
        .allow_auto_conversion("RecordA")
        .close();

    static Type::Record descA_rec = Type::Record( "RecordA", "desc.")
        .derive_from(a_rec)
		.copy_keys(inrec)
	    .close();

    static Type::Record root_rec_with_abstract = Type::Record( "RootRecordWithAbstract", "desc.")
        .declare_key("format", a_rec, Type::Default::obligatory(), "")
	    .close();

    { // JSON format, read record
        stringstream ss("{ format: {} }");
        read_stream(ss, root_rec);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecord", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ("InRecord", storage_->get_item(1)->get_item(0)->get_string() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // YAML format, read record
        stringstream ss(input_yaml_empty_rec_without_tag);
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecord", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ("InRecord", storage_->get_item(1)->get_item(0)->get_string() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // JSON format, read abstract
        stringstream ss("{ format: { TYPE=\"RecordA\" } }");
        read_stream(ss, root_rec_with_abstract);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecordWithAbstract", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "RecordA");
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // YAML format, read abstract
        stringstream ss(input_yaml_empty_rec_with_tag);
        read_stream(ss, root_rec_with_abstract, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecordWithAbstract", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "RecordA");
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // YAML format, read abstract with autoconversion
        stringstream ss(input_yaml_empty_rec_without_tag);
        read_stream(ss, root_rec_with_abstract, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecordWithAbstract", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "RecordA");
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // YAML format, read record with empty keys
        stringstream ss(input_yaml_rec_with_empty_scalar);
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ("RootRecord", storage_->get_item(0)->get_string());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ("InRecord", storage_->get_item(1)->get_item(0)->get_string() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(1)->is_null() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

}



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
    	.declare_key("sel_key", sel_type, Type::Default("\"two\""),"")
    	.declare_key("double_key", Type::Double(), Type::Default("1.23"),"")
    	.declare_key("str_key", Type::String(), Type::Default("\"ahoj\""),"")
    	.declare_key("array_key", Type::Array( Type::Integer() ), Type::Default("123"), "")
	    .declare_key("rec_key", sub_rec, Type::Default("321"), "")
		.close();

    rec_type.finish();
    sub_rec.finish();

    {
        stringstream ss("{ }");
        read_stream(ss, rec_type);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(8, storage_->get_array_size());
        EXPECT_EQ("SomeRec", storage_->get_item(0)->get_string());
        EXPECT_EQ(4, storage_->get_item(1)->get_int() );
        EXPECT_TRUE( storage_->get_item(2)->get_bool() );
        EXPECT_EQ(2, storage_->get_item(3)->get_int() );
        EXPECT_EQ(1.23, storage_->get_item(4)->get_double() );
        EXPECT_EQ("ahoj", storage_->get_item(5)->get_string() );
        EXPECT_EQ(123 , storage_->get_item(6)->get_item(0)->get_int() );
        EXPECT_EQ("SubRecord", storage_->get_item(7)->get_item(0)->get_string() );
        EXPECT_FALSE( storage_->get_item(7)->get_item(1)->get_bool() );
        EXPECT_EQ(321 , storage_->get_item(7)->get_item(2)->get_int() );
    }
}


TEST_F(InputReaderToStorageTest, storage_transpose) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Selection sel_type = Type::Selection("IntSelection")
    	.add_value(10,"ten","")
    	.add_value(1,"one","")
    	.add_value(2,"two","")
		.close();

    static Type::Record subrec_2 = Type::Record("TransposeSubrec2","")
		.declare_key("key_a", Type::String(), "")
		.declare_key("key_b", Type::Bool(), "")
		.close();

    static Type::Record subrec_1 = Type::Record("TransposeSubrec1","")
		.declare_key("one", Type::Integer(), "")
		.declare_key("two", Type::Array( Type::Integer(),1,10 ), "")
		.declare_key("three", subrec_2, "")
		.declare_key("four", Type::Double(), "")
		.declare_key("five", sel_type, "")
		.close();

    static Type::Record root_record = Type::Record("RootTransposeRec","")
		.declare_key("set", Type::Array(subrec_1,1,3), "")
		.declare_key("default", Type::Bool(), "")
		.close();

    {   // Try one correct type
        stringstream ss("{ set={ one=[1,2,3], two=[2,3], three={key_a=[\"A\",\"B\",\"C\"], key_b=[false,true,false]}, "
        		"four=[1.5, 2.5, 3.5], five=[\"one\",\"two\",\"ten\"] }, default=true }");
        read_stream(ss, root_record);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(3,   storage_->get_item(1)->get_array_size());
        EXPECT_TRUE(   storage_->get_item(2)->get_bool());
        StorageBase *substorage = storage_->get_item(1)->get_item(0);
        EXPECT_EQ(6,   substorage->get_array_size());
        EXPECT_EQ(1,   substorage->get_item(1)->get_int());
        EXPECT_EQ(2,   substorage->get_item(2)->get_array_size());
        EXPECT_EQ(3,   substorage->get_item(3)->get_array_size());
        EXPECT_EQ("A", substorage->get_item(3)->get_item(1)->get_string());
        EXPECT_FALSE(  substorage->get_item(3)->get_item(2)->get_bool());
        EXPECT_EQ(1.5, substorage->get_item(4)->get_double());
        EXPECT_EQ(1,   substorage->get_item(5)->get_int());
        substorage = storage_->get_item(1)->get_item(1);
        EXPECT_EQ(6,   substorage->get_array_size());
        EXPECT_EQ(2,   substorage->get_item(1)->get_int());
        EXPECT_EQ("B", substorage->get_item(3)->get_item(1)->get_string());
        EXPECT_TRUE(   substorage->get_item(3)->get_item(2)->get_bool());
        EXPECT_EQ(2.5, substorage->get_item(4)->get_double());
        EXPECT_EQ(2,   substorage->get_item(5)->get_int());
        substorage = storage_->get_item(1)->get_item(2);
        EXPECT_EQ(6,   substorage->get_array_size());
        EXPECT_EQ(3,   substorage->get_item(1)->get_int());
        EXPECT_EQ("C", substorage->get_item(3)->get_item(1)->get_string());
        EXPECT_FALSE(  substorage->get_item(3)->get_item(2)->get_bool());
        EXPECT_EQ(3.5, substorage->get_item(4)->get_double());
        EXPECT_EQ(10,  substorage->get_item(5)->get_int());
    }

    {   // Try other correct type - automatic conversion to array with one element
        stringstream ss("{ set={ one=1, two=[2,3], three={key_a=\"A\", key_b=false}, four=1.5, five=\"one\" }, default=true }");
        read_stream(ss, root_record);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(1, storage_->get_item(1)->get_array_size());
        StorageBase *substorage = storage_->get_item(1)->get_item(0);
        EXPECT_EQ(6,   substorage->get_array_size());
        EXPECT_EQ(1,   substorage->get_item(1)->get_int());
        EXPECT_EQ(2,   substorage->get_item(2)->get_array_size());
        EXPECT_EQ(3,   substorage->get_item(3)->get_array_size());
        EXPECT_EQ("A", substorage->get_item(3)->get_item(1)->get_string());
        EXPECT_FALSE(  substorage->get_item(3)->get_item(2)->get_bool());
        EXPECT_EQ(1.5, substorage->get_item(4)->get_double());
    }

    {   // Incorrect type - empty sub-array
    	stringstream ss("{ set={ one=[], two=[2,3], three={key_a=[\"A\",\"B\",\"C\"], key_b=false}, four=[1.5, 2.5, 3.5], five=\"one\" }, default=true }");
        EXPECT_THROW_WHAT( {read_stream(ss, root_record);}, ReaderInternalBase::ExcInputError, "Empty array during transpose auto-conversion.");
    }

    {   // Incorrect type - unequal sizes of sub-arrays
    	stringstream ss("{ set={ one=[1,2,3], two=[2,3], three={key_a=[\"A\",\"B\"], key_b=false}, four=[1.5, 2.5, 3.5], five=\"one\" }, default=true }");
        EXPECT_THROW_WHAT( {read_stream(ss, root_record);}, ReaderInternalBase::ExcInputError, "Unequal sizes of sub-arrays during transpose auto-conversion");
    }

    {   // Incorrect type - transposition of Record is not allowed
        stringstream ss("{ set={ one=[1,2], two=[2,3], three=[{key_a=\"A\", key_b=false}, {key_a=\"B\", key_b=false}], four=1.5, five=\"one\" }, default=true }");
        EXPECT_THROW_WHAT( {read_stream(ss, root_record);}, ReaderInternalBase::ExcInputError, "The value should be 'JSON object', but we found:.* 'JSON array'");
    }

    {   // Incorrect type - array is greater than size limit
        stringstream ss("{ set={ one=[1,2,3,4], two=[2,3], three={key_a=\"B\", key_b=false}, four=1.5, five=\"one\" }, default=true }");
        EXPECT_THROW_WHAT( {read_stream(ss, root_record);}, ReaderInternalBase::ExcInputError, "Result of transpose auto-conversion do not fit the size 4 of the Array");
    }

}


const string input_yaml_noautoconversion = R"YAML(
pressure: !DescendantB
  value: 0.5
)YAML";

const string input_json_noautoconversion = R"JSON(
{ 
  pressure: {
    TYPE = "DescendantB",
    value: 0.5
  }
}
)JSON";

TEST_F(InputReaderToStorageTest, Abstract_auto_conversion) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    static Type::Abstract abstract = Type::Abstract("AbstractOfRec", "Base of records.")
        .allow_auto_conversion("DescendantA")
        .close();

    static Type::Record desc_a = Type::Record( "DescendantA", "descendant A")
        .derive_from(abstract)
		.declare_key("value", Type::Double(), Type::Default::obligatory(), "value")
		.allow_auto_conversion("value")
	    .close();

    static Type::Record desc_b = Type::Record( "DescendantB", "descendant B")
        .derive_from(abstract)
		.declare_key("value", Type::Double(), Type::Default::obligatory(), "value")
		.declare_key("b_val", Type::Integer(), Type::Default::optional(), "value")
		.allow_auto_conversion("value")
	    .close();

    static Type::Record desc_c = Type::Record( "DescendantC", "descendant C")
        .derive_from(abstract)
		.declare_key("value", Type::Double(), Type::Default::obligatory(), "value")
		.declare_key("c_val", Type::String(), Type::Default::optional(), "value")
	    .close();

    static Type::Record root_rec = Type::Record("RootOfAutoconversionTest", "Root of IST.")
        .declare_key("pressure", abstract, Type::Default::obligatory(), "presure")
        .close();

    { // YAML format, DescendantA with autoconversion of Abstract and Record
        stringstream ss("pressure: 0.5");
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(2, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantA");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
    }

    { // JSON format, DescendantA with autoconversion of Abstract and Record
        stringstream ss("{ pressure = 0.5 }");
        read_stream(ss, root_rec, FileFormat::format_JSON);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(2, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantA");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
    }

    { // YAML format, DescendantA only with autoconversion of Record
        stringstream ss("pressure: !DescendantA 0.5");
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(2, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantA");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
    }

    { // YAML format, DescendantB with autoconversion
        stringstream ss("pressure: !DescendantB 0.5");
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantB");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // YAML format, DescendantC with autoconversion
        stringstream ss("pressure: !DescendantC 0.5");
        EXPECT_THROW_WHAT( {read_stream(ss, root_rec, FileFormat::format_YAML);}, ReaderInternalBase::ExcInputError,
        		"The value should be 'YAML map', but we found");
    }

    { // YAML format, full output
        stringstream ss(input_yaml_noautoconversion);
        read_stream(ss, root_rec, FileFormat::format_YAML);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantB");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

    { // JSON format, full output
        stringstream ss(input_json_noautoconversion);
        read_stream(ss, root_rec, FileFormat::format_JSON);

        EXPECT_NE((void *)NULL, storage_);
        EXPECT_EQ(2, storage_->get_array_size());
        EXPECT_EQ(3, storage_->get_item(1)->get_array_size());
        EXPECT_EQ(storage_->get_item(1)->get_item(0)->get_string(), "DescendantB");
        EXPECT_DOUBLE_EQ(0.5, storage_->get_item(1)->get_item(1)->get_double() );
        EXPECT_TRUE(storage_->get_item(1)->get_item(2)->is_null() );
    }

}


