/*
 * input_type_test.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 */


#include <gtest/gtest.h>

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
 * Test Array class.
 */
TEST(InputTypeArray, all_methods) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    // construction
    Array arr_int(Integer(),1,8);
    Array arr_arr_dbl( Array( Double() ));

    Record rec_2("record_type_2", "desc");
    rec_2.finish();

    Array arr_rec_shared_ptr( rec_2 );

    Selection<int> sel("Singular set.");
    sel.finish();

    Array arr_of_sel( sel );

    // get_sub_type
    EXPECT_EQ( rec_2, arr_rec_shared_ptr.get_sub_type());

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

}

/**
 * Test all simple scalar types.
 */
TEST(InputTypeScalar, all_types) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    EXPECT_EQ( FileName(output_file), FileName(output_file) );
    EXPECT_NE( FileName(output_file), FileName(input_file) );
    EXPECT_NE( FileName(output_file), Integer() );
    EXPECT_EQ( output_file, FileName(output_file).get_file_type() );
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

    Selection<enum Colors> *sel1= new Selection<enum Colors>("Colors");
    Selection<enum Colors> sel2;
    sel2=*sel1;


    sel1->add_value(blue, "blue");
    sel2.add_value(white,"white","White color");
    sel1->add_value(black,"black");
    sel2.add_value(red,"red");
    EXPECT_DEATH( {sel1->add_value(green,"red");}, "already exists in Selection:");
    EXPECT_DEATH( {sel2.add_value(green,"red");}, "already exists in Selection:");
    EXPECT_DEATH( {sel1->add_value(blue,"blue1");}, "conflicts with value");
    EXPECT_DEATH( {sel2.add_value(blue,"blue1");}, "conflicts with value");
    delete sel1;

    sel2.add_value(green,"green");
    sel2.finish();
    EXPECT_DEATH( {sel2.add_value(yellow,"y");}, "in finished Selection type:");

    Selection<int> sel3;
    EXPECT_DEATH( {sel3.add_value(1,"one");}, "to empty selection handle." );
    // getter methods
    EXPECT_TRUE( sel2.has_name("blue") );
    EXPECT_FALSE( sel2.has_name("xblue") );
    EXPECT_TRUE( sel2.has_value(blue) );
    EXPECT_TRUE( sel2.has_value(black) );
    EXPECT_FALSE( sel2.has_value(yellow) );

    EXPECT_EQ( 45, sel2.name_to_int("black") );
    EXPECT_THROW( {sel2.name_to_int("xblack");}, SelectionBase::ExcSelectionKeyNotFound );


    // Integer selection
    Selection<int> int_sel("Integer selection");
    int_sel.add_value(10, "ten");
    int_sel.add_value(0,"zero");
    int_sel.finish();

    EXPECT_EQ(10, int_sel.name_to_int("ten"));
}


/**
 * Test documentation output.
 */
TEST(InputTypeDocumentation, whole_tree) {
using namespace Input::Type;

Record output_record("OutputRecord",
        "Information about one file for field data.");
{
    output_record.declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::optional),
            "File for output stream.");

    output_record.declare_key("digits",Integer(0,8), DefaultValue("8"),
            "Number of digits used for output double values into text output files.");
    output_record.declare_key("compression", Bool(),
            "Whether to use compression of output file.");

    output_record.declare_key("start_time", Double(0.0),
            "Simulation time of first output.");
    output_record.declare_key("data_description", String(),DefaultValue(),
            "");
    output_record.finish();
} // delete local variables

Record array_record("RecordOfArrays",
         "Long description of record.\n"
         "Description could have more lines"
         );
{
 Array array_of_int(Integer(0), 5, 100 );

 array_record.declare_key("array_of_5_ints", array_of_int,
         "Some bizare array.");
 array_record.declare_key("array_of_str", Array( String() ), DefaultValue(),
         "Desc. of array");
 array_record.declare_key("array_of_str_1", Array( String() ), DefaultValue(),
             "Desc. of array");
 array_record.finish();
}


 Record record_record("RecordOfRecords",
         "Long description of record.\n"
         "Description could have more lines"
         );


 {
     Record other_record("OtherRecord","desc");
     other_record.finish();

     record_record.declare_key("sub_rec_1", other_record, "key desc");

     // recursion
     //record_record->declare_key("sub_rec_2", record_record, "Recursive key.");

     record_record.finish();
 }

 Selection<enum Colors> sel("Colors");
 {
     sel.add_value(blue, "blue");
     sel.add_value(white,"white","White color");
     sel.add_value(black,"black");
     sel.add_value(red,"red");
     sel.add_value(green,"green");
     sel.finish();
 }

 Record main("MainRecord", "The main record of flow.");
 main.declare_key("array_of_records", Array(output_record), "Array of output streams.");
 main.declare_key("record_record", record_record, "no comment on record_record");
 main.declare_key("color", sel, "My favourite color.");
 main.declare_key("color1", sel, "My second favourite color.");
 main.declare_key("array_record", array_record, "no commment on array_record");
 main.finish();


 main.documentation(cout, true);
}
