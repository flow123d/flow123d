/*
 * input_type_test.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <input/input_type.hh>



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


    // from_default methods
    // Bool
    EXPECT_TRUE( Bool().from_default("true") );
    EXPECT_FALSE( Bool().from_default("false") );
    EXPECT_THROW_WHAT( { Bool().from_default("yes"); }, ExcWrongDefault,
            "Default value 'yes' do not match type: 'Bool';" );

    // Integer
    EXPECT_EQ(10, Integer().from_default("10") );
    EXPECT_THROW_WHAT( { Integer(0,4).from_default("10"); }, ExcWrongDefault,
            "Default value '10' do not match type: 'Integer';" );
    EXPECT_THROW_WHAT( { Integer().from_default("yes"); }, ExcWrongDefault,
            "Default value 'yes' do not match type: 'Integer';" );

    // Double
    EXPECT_EQ(3.14, Double().from_default("3.14") );
    EXPECT_EQ(-5.67e-23, Double().from_default("-5.67E-23") );
    EXPECT_THROW_WHAT( { Double(0,4.4).from_default("-1e-10"); }, ExcWrongDefault,
            "Default value .* do not match type: 'Double';" );
    EXPECT_THROW_WHAT( { Double().from_default("-3.6t5"); }, ExcWrongDefault,
            "Default value .* do not match type: 'Double';" );


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
    Array arr_int(Integer(),1,8);
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
    EXPECT_NE( arr_int, Array( Double() ) );
    EXPECT_EQ( arr_int, Array( Integer() ) );
    EXPECT_NE( arr_int, rec_2 );

    // match_size
    EXPECT_TRUE( arr_int.match_size(1) );
    EXPECT_TRUE( arr_int.match_size(3) );
    EXPECT_TRUE( arr_int.match_size(8) );
    EXPECT_FALSE( arr_int.match_size(9) );

    // type_name()
    EXPECT_EQ( arr_int.type_name(), Array(Integer()).type_name() );

    // valid default
    // no death
    arr_int.valid_default("100");
    EXPECT_THROW_WHAT( {arr_int.valid_default("1.5");}, ExcWrongDefault,
            "Default value '1.5' do not match type:" );
    EXPECT_THROW_WHAT( {Array( Double(), 2 ).valid_default("3.2"); }, ExcWrongDefault,
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
    Selection int_sel("Integer selection");
    int_sel.add_value(10, "ten");
    int_sel.add_value(0,"zero");
    int_sel.close();


    // Selection defaults
    EXPECT_EQ(10, int_sel.from_default("ten") );
    EXPECT_EQ(0, int_sel.from_default("zero") );
    EXPECT_THROW_WHAT( { int_sel.from_default("two"); }, ExcWrongDefault,
            "Default value .* do not match type: " );
    EXPECT_THROW_WHAT( { int_sel.from_default("10"); }, ExcWrongDefault,
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
