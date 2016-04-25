/*
 * input_type_test.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <input/input_type.hh>
#include <input/attribute_lib.hh>
#include <input/type_attribute_lib.hh>



/**
 * todo:
 * - testovat ze je od kazdeho rRecordu jen jedna instance
 * - testovat  ze Scalarni typy kontroluji platnost defaultnich hodnot
 */

/**
 * Check is boost type traits can distinguish int and enum.
 */
#include <boost/type_traits.hpp>

enum x_enum {
    one=1,
    two=2
};
typedef x_enum XX;

TEST(BoostTypeTraits, Enum) {
    EXPECT_TRUE( boost::is_enum<XX>::value);
    EXPECT_FALSE( boost::is_integral<XX>::value);
};


class PublicTypeBase : public Input::Type::TypeBase {
public:
    using Input::Type::TypeBase::key_hash;
    using Input::Type::TypeBase::is_valid_identifier;
};

/**
 * Test TypeBase class
 */
TEST(InputTypeTypeBase, static_methods) {
    using namespace Input::Type;
    // can call static methods
    EXPECT_NE( PublicTypeBase::key_hash("Ahoj"), PublicTypeBase::key_hash("Cau") );
    EXPECT_TRUE( PublicTypeBase::is_valid_identifier("00ahoj_0123456789"));
    EXPECT_FALSE( PublicTypeBase::is_valid_identifier("Ahoj"));
    EXPECT_FALSE( PublicTypeBase::is_valid_identifier("$%"));
    EXPECT_FALSE( PublicTypeBase::is_valid_identifier("a h o j"));
}


/**
 * Test all simple scalar types.
 */
TEST(InputTypeScalar, all_types) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";
	using namespace Input::Type;


    // Default::check_validity methods
    // Bool
    EXPECT_TRUE( Default("true").check_validity( boost::make_shared<Bool>() ) );
    EXPECT_TRUE( Default("false").check_validity( boost::make_shared<Bool>() ) );
    EXPECT_FALSE( Default::obligatory().check_validity( boost::make_shared<Bool>()) );
    EXPECT_THROW_WHAT( { Default("yes").check_validity( boost::make_shared<Bool>() ); }, ExcWrongDefault,
            "Default value 'yes' do not match type: 'Bool';" );

    // Integer
    EXPECT_TRUE( Default("10").check_validity( boost::make_shared<Integer>() ) );
    EXPECT_THROW_WHAT( { Default("10").check_validity( boost::make_shared<Integer>(0,4) ); }, ExcWrongDefault,
            "Default value '10' do not match type: 'Integer';" );
    EXPECT_THROW_WHAT( { Default("yes").check_validity( boost::make_shared<Integer>() ); }, ExcWrongDefault,
            "Default value 'yes' do not match type: 'Integer';" );

    // Double
    EXPECT_TRUE( Default("3.14").check_validity( boost::make_shared<Double>() ) );
    EXPECT_TRUE( Default("-5.67E-23").check_validity( boost::make_shared<Double>() ) );
    EXPECT_THROW_WHAT( { Default("-1e-10").check_validity( boost::make_shared<Double>(0,4.4) ); }, ExcWrongDefault,
            "Default value .* do not match type: 'Double';" );
    EXPECT_THROW_WHAT( { Default("-3.6t5").check_validity( boost::make_shared<Double>() ); }, ExcWrongDefault,
            "Default value .* do not match type: 'Double';" );

    // String
    EXPECT_TRUE( Default("\"ahoj\"").check_validity( boost::make_shared<String>() ) );
    EXPECT_THROW_WHAT( { Default("ahoj").check_validity( boost::make_shared<String>() ); }, ExcWrongDefault,
            "Default value .* do not match type: 'String';" );


    // test equivalence operator
    EXPECT_EQ( FileName::output(), FileName::output() );
    EXPECT_NE( FileName::output(), FileName::input() );
    EXPECT_NE( FileName::output(), Integer() );

    // test getter for file type
    EXPECT_EQ( FilePath::output_file, FileName::output().get_file_type() );


}

/**
 * Test Array class.
 */
TEST(InputTypeArray, all_methods) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // construction
    boost::shared_ptr<Array> arr_int = boost::make_shared<Array>(Integer(), 1, 8);
    Array arr_arr_dbl( Array( Double() ));

    Record rec_2("record_type_2", "desc");
    rec_2.close();

    Array arr_rec_shared_ptr( rec_2 );

    Selection sel("Singular set.");
    sel.close();

    Array arr_of_sel( sel );

    Input::Type::TypeBase::lazy_finish();

    // get_sub_type
    EXPECT_EQ( rec_2, arr_rec_shared_ptr.get_sub_type()); // boost::smart_ptr assert fails

    // operator ==
    EXPECT_NE( *arr_int, Array( Double() ) );
    EXPECT_EQ( *arr_int, Array( Integer() ) );
    EXPECT_NE( *arr_int, rec_2 );

    // match_size
    EXPECT_TRUE( arr_int->match_size(1) );
    EXPECT_TRUE( arr_int->match_size(3) );
    EXPECT_TRUE( arr_int->match_size(8) );
    EXPECT_FALSE( arr_int->match_size(9) );

    // type_name()
    EXPECT_EQ( arr_int->type_name(), Array(Integer()).type_name() );

    // valid default
    // no death
    Default("100").check_validity( arr_int );
    EXPECT_THROW_WHAT( {Default("1.5").check_validity( arr_int );}, ExcWrongDefault,
            "Default value '1.5' do not match type:" );
    boost::shared_ptr<Array> arr_double = boost::make_shared<Array>(Double(), 2);
    EXPECT_THROW_WHAT( { Default("3.2").check_validity( arr_double ); }, ExcWrongDefault,
                  "Default value '3.2' do not match type: 'array_of_Double';"
                 );

}




/**
 * Test Selection class.
 */
enum Colors {
    blue,
    white=300,
    black=45,
    red,
    green,
    yellow
};

TEST(InputTypeSelection, construction) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Selection *sel1= new Selection("Colors");
    Selection sel2=*sel1;

    sel1->add_value(blue, "blue");
    sel2.add_value(white,"white","White color");
    sel1->add_value(black,"black");
    sel2.add_value(red,"red");
    EXPECT_THROW_WHAT( {sel1->add_value(green,"red");}, ExcXprintfMsg, "already exists in Selection:");
    EXPECT_THROW_WHAT( {sel2.add_value(green,"red");}, ExcXprintfMsg, "already exists in Selection:");
    EXPECT_THROW_WHAT( {sel1->add_value(blue,"blue1");}, ExcXprintfMsg, "conflicts with value");
    EXPECT_THROW_WHAT( {sel2.add_value(blue,"blue1");}, ExcXprintfMsg, "conflicts with value");


    sel2.add_value(green,"green");
    sel2.close();
    EXPECT_THROW_WHAT( {sel2.add_value(yellow,"y");}, ExcXprintfMsg, "in finished Selection type:");

    Selection sel3;
    EXPECT_TRUE( sel3.is_finished());
    EXPECT_EQ("EmptySelection", sel3.type_name());
    EXPECT_THROW_WHAT( {sel3.add_value(1,"one");}, ExcXprintfMsg, "in finished Selection type:");
    // getter methods
    EXPECT_TRUE( sel2.has_name("blue") );
    EXPECT_FALSE( sel2.has_name("xblue") );
    EXPECT_TRUE( sel2.has_value(blue) );
    EXPECT_TRUE( sel2.has_value(black) );
    EXPECT_FALSE( sel2.has_value(yellow) );

    EXPECT_EQ( 45, sel2.name_to_int("black") );
    EXPECT_THROW( {sel2.name_to_int("xblack");}, Selection::ExcSelectionKeyNotFound );


    // Integer selection
    Selection int_sel = Selection("Integer selection")
    	.add_value(10, "ten")
    	.add_value(0,"zero")
    	.close();
    boost::shared_ptr<Selection> sel_ptr = boost::make_shared<Selection>(int_sel);


    // Selection defaults
    EXPECT_TRUE( Default("\"ten\"").check_validity( sel_ptr ) );
    EXPECT_TRUE( Default("\"zero\"").check_validity( sel_ptr ) );
    EXPECT_THROW_WHAT( { Default("\"two\"").check_validity( sel_ptr ); }, ExcWrongDefault,
            "Default value .* do not match type: " );
    EXPECT_THROW_WHAT( { Default("10").check_validity( sel_ptr ); }, ExcWrongDefault,
            "Default value .* do not match type: " );

    EXPECT_EQ(10, int_sel.name_to_int("ten"));

    delete sel1;
}


/**
 * Test documentation output.
 */
TEST(InputTypeDocumentation, whole_tree) {
using namespace Input::Type;

Record output_record("OutputRecord",
        "Information about one file for field data.");
{
    output_record.declare_key("file", FileName::output(), Default::optional(),
            "File for output stream.");

    output_record.declare_key("digits",Integer(0,8), Default("8"),
            "Number of digits used for output double values into text output files.");
    output_record.declare_key("compression", Bool(),
            "Whether to use compression of output file.");

    output_record.declare_key("start_time", Double(0.0),
            "Simulation time of first output.");
    output_record.declare_key("data_description", String(), Default::optional(),
            "");
    output_record.close();
} // delete local variables

Record array_record("RecordOfArrays",
         "Long description of record.\n"
         "Description could have more lines"
         );
{
 Array array_of_int(Integer(0), 5, 100 );

 array_record.declare_key("array_of_5_ints", array_of_int,
         "Some bizare array.");
 array_record.declare_key("array_of_str", Array( String() ), Default::optional(),
         "Desc. of array");
 array_record.declare_key("array_of_str_1", Array( String() ), Default::optional(),
             "Desc. of array");
 array_record.close();
}


 Record record_record("RecordOfRecords",
         "Long description of record.\n"
         "Description could have more lines"
         );


 {
     Record other_record("OtherRecord","desc");
     other_record.close();

     record_record.declare_key("sub_rec_1", other_record, "key desc");

     // recursion
     //record_record->declare_key("sub_rec_2", record_record, "Recursive key.");

     record_record.close();
 }

 Selection sel("Colors");
 {
     sel.add_value(blue, "blue");
     sel.add_value(white,"white","White color");
     sel.add_value(black,"black");
     sel.add_value(red,"red");
     sel.add_value(green,"green");
     sel.close();
 }

 Record main("MainRecord", "The main record of flow.");
 main.declare_key("array_of_records", Array(output_record), "Array of output streams.");
 main.declare_key("record_record", record_record, "no comment on record_record");
 main.declare_key("color", sel, "My favourite color.");
 main.declare_key("color1", sel, "My second favourite color.");
 main.declare_key("array_record", array_record, "no commment on array_record");
 main.close();

// main.documentation(cout, TypeBase::full_after_record);

}


class InputTypeAttributesTest : public testing::Test, public Input::Type::Record {
public:
	InputTypeAttributesTest() : Input::Type::Record("rec", "Some record.") {}
protected:

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };

};

TEST_F(InputTypeAttributesTest, base_test) {
	std::map<std::string, json_string>::iterator it;

	this->add_attribute("attr_1", "\"some attribute\"");
	this->add_attribute("attr_2", "\"other attribute\"");
	this->add_attribute("numeric", "\"10\"");
	this->add_attribute("pair", "[\"0\", \"50\"]");
	this->add_attribute("float_point", "\"0.5\"");
	this->add_attribute(FlowAttributes::parameters(), "[\"a\", \"b\", \"c\"]");
	this->add_attribute(Input::Type::Attributes::obsolete(), "\"true\"");

	EXPECT_TRUE( (it=attributes_->find("attr_1")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "\"some attribute\"" );
	EXPECT_TRUE( (it=attributes_->find("numeric")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "\"10\"" );
	EXPECT_TRUE( (it=attributes_->find("pair")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "[\"0\", \"50\"]" );
	EXPECT_TRUE( (it=attributes_->find("float_point")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "\"0.5\"" );
	EXPECT_TRUE( (it=attributes_->find("parameters")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "[\"a\", \"b\", \"c\"]" );
	EXPECT_TRUE( (it=attributes_->find("obsolete")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "\"true\"" );

	this->add_attribute("numeric", "\"5\"");
	EXPECT_TRUE( (it=attributes_->find("numeric")) != attributes_->end() );
	EXPECT_STREQ( it->second.c_str(), "\"5\"" );

    EXPECT_THROW_WHAT( { this->add_attribute("invalid_attr", "non quotation attribute"); }, ExcXprintfMsg,
            "Invalid JSON format of attribute 'invalid_attr'." );
}

/*TEST(InputTypeAttributes, complete_test) {
	using namespace Input::Type;

	static Selection sel = Selection("Bc_types", "Boundary conditions.")
		.add_value(0, "none")
		.add_value(1, "Newton")
		.add_value(2, "Neumann")
		.add_value(3, "Dirichlet")
		.add_value(4, "Robin")
		.close();

	static Record rec1 = Record("rec_1", "Time step record.")
		.declare_key("time_step", Double(), Default::obligatory(), "Time step value.")
		.close();

	static Record rec2 = Record("rec_2", "Step record.")
		.declare_key("step_iter", Integer(1, 100), Default::obligatory(), "Count of steps.")
		.close();

	static Abstract a_rec = Abstract("a_rec", "Step Abstract")
		.close();
	a_rec.add_child(rec1);
	a_rec.add_child(rec2);

	static Record main_rec = Record("main_rec", "Main Record.")
		.declare_key("bc", sel, Default::obligatory(), "Selection of boundary conditions.")
		.declare_key("idx", Array( Integer() ), Default::obligatory(), "List of indexes.")
		.declare_key("time", Double(), Default::obligatory(), "Start time value.")
		.declare_key("step", a_rec, Default::obligatory(), "Start time value.")
		.declare_key("desc", String(), Default::optional(), "Description of problem.")
		.declare_key("file", FileName::output(), Default::optional(), "File for output stream.")
		.close();

	TypeBase::lazy_finish();

	{
		std::stringstream ss;
		main_rec.write_attributes(ss);
		string str = ss.str();
		//EXPECT_TRUE( str.find("\"full_name\" : \"main_rec\"") != std::string::npos );
	}
	cout << "---------------" << endl;
	main_rec.write_attributes(cout);
	cout << "---------------" << endl;
	a_rec.write_attributes(cout);
	cout << "---------------" << endl;
	sel.write_attributes(cout);
	cout << "---------------" << endl;
}*/
