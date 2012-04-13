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


/**
 * Test Array class.
 */
TEST(InputTypeArray, construction) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Array arr_int(Integer());
    Array arr_arr_dbl( Array( Double() ));

    Record rec_1("record_type_1","desc");
    boost::shared_ptr<Record> rec_2 = boost::make_shared<Record>("record_type_2", "desc");


    EXPECT_DEATH( {Array arr_rec_ref( Record("subrec_type", "desc") ); },
                 "Complex type .* shared_ptr."
                );

    Array arr_rec_shared_ptr( rec_2 );
}

/**
 * Test Selection class.
 */
enum Colors {
    blue,
    white=300,
    black=45,
    red,
    green
};

TEST(InputTypeSelection, construction) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Selection<enum Colors> sel("Colors");
    sel.add_value(blue, "blue");
    sel.add_value(white,"white","White color");
    sel.add_value(black,"black");
    sel.add_value(red,"red");
    EXPECT_DEATH( {sel.add_value(green,"red");}, "Existing name in declaration");
    EXPECT_DEATH( {sel.add_value(blue,"blue1");}, "Existing value in declaration");
    sel.add_value(green,"green");

}

/**
 * Test Record class.
 */
TEST(InputTypeRecord, declare_key_scalars) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FileName
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   static Record rec("SomeRecord", "desc.");

   rec.declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::read_time), "desc.");

   boost::shared_ptr<Integer> digits_type=boost::make_shared<Integer>((int)0, (int)8);
   rec.declare_key("digits",digits_type, DefaultValue("8"), "desc.");

   rec.declare_key("compression", Bool(),"desc.");

   boost::shared_ptr<Double> time_type=boost::make_shared<Double>(0.0);
   rec.declare_key("start_time", time_type,"desc.");

   rec.declare_key("data_description", String(),DefaultValue(),"");

   EXPECT_DEATH( { rec.declare_key("data_description", String(),"");},
                "Redeclaration of key:");

   EXPECT_DEATH( { DefaultValue d(DefaultValue::declaration);},
                "Can not construct DefaultValue with type 'declaration' without providing the default value.");

   enum Colors {
       white, black, red
   };

   MAKE_Input_Type_Selection(Colors, sel)("Color selection");
   sel->add_value(black, "black");
   sel->add_value(red, "red");

   rec.declare_key("plot_color", sel, "Color to plot the fields in file.");

   EXPECT_DEATH( { rec.declare_key("x", *sel, "desc.");},
                "Complex type .* shared_ptr.");

   // should fail at compile time
   // output_record.declare_key("xx", 10, "desc");

   //boost::shared_ptr<int> int_ptr = boost::make_shared<int>(8);
   //output_record.declare_key("xx", int_ptr, "desc");
}

TEST(InputTypeRecord, declare_key_arrays) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


   static Record array_record("RecordOfArrays", "desc.");

   // array type passed through shared_ptr
    boost::shared_ptr<Array> array_of_int = boost::make_shared<Array>(Integer(0), 5, 100 );
    array_record.declare_key("array_of_5_ints", array_of_int,"Some bizare array.");

    // array type passed by reference
    array_record.declare_key("array_of_str", Array( String() ), DefaultValue(),"Desc. of array");
    array_record.declare_key("array_of_str_1", Array( String() ), DefaultValue(), "Desc. of array");


    ASSERT_DEATH( {array_record.declare_key("some_key", Array( String() ), DefaultValue("10"), ""); },
                  "Can not provide default value for non scalar type."
                 );

}

TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


    boost::shared_ptr<Record> record_record=boost::make_shared<Record>("RecordOfRecords", "");

    ASSERT_DEATH( {record_record->declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
                  "Complex type .* shared_ptr."
                  );

    boost::shared_ptr<Record> other_record=boost::make_shared<Record>("OtherRecord","desc");
    record_record->declare_key("sub_rec_1", other_record, "key desc");

    // recursion
    record_record->declare_key("sub_rec_2", record_record, "desc.");

}

/**
 * Test documentation output.
 */
TEST(InputTypeDocumentation, whole_tree) {
using namespace Input::Type;

boost::shared_ptr<Record> output_record=boost::make_shared<Record>("OutputRecord",
        "Information about one file for field data.");
{
    output_record->declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::read_time),
            "File for output stream.");
    boost::shared_ptr<Integer> digits_type=boost::make_shared<Integer>((int)0, (int)8);
    output_record->declare_key("digits",digits_type, DefaultValue("8"),
            "Number of digits used for output double values into text output files.");
    output_record->declare_key("compression", Bool(),
            "Whether to use compression of output file.");
    boost::shared_ptr<Double> time_type=boost::make_shared<Double>(0.0);
    output_record->declare_key("start_time", time_type,
            "Simulation time of first output.");
    output_record->declare_key("data_description", String(),DefaultValue(),
            "");
} // delete local variables

boost::shared_ptr<Record> array_record=boost::make_shared<Record>("RecordOfArrays",
         "Long description of record.\n"
         "Description could have more lines"
         );
{
 boost::shared_ptr<Array> array_of_int = boost::make_shared<Array>(Integer(0), 5, 100 );
 array_record->declare_key("array_of_5_ints", array_of_int,
         "Some bizare array.");
 array_record->declare_key("array_of_str", Array( String() ), DefaultValue(),
         "Desc. of array");
 array_record->declare_key("array_of_str_1", Array( String() ), DefaultValue(),
             "Desc. of array");
}


 boost::shared_ptr<Record> record_record= boost::make_shared<Record>("RecordOfRecords",
         "Long description of record.\n"
         "Description could have more lines"
         );


 {
     boost::shared_ptr<Record> other_record=boost::make_shared<Record>("OtherRecord","desc");
     record_record->declare_key("sub_rec_1", other_record, "key desc");

     // recursion
     record_record->declare_key("sub_rec_2", record_record, "Recursive key.");
 }

 MAKE_Input_Type_Selection(enum Colors, sel)("Colors");
 {
     sel->add_value(blue, "blue");
     sel->add_value(white,"white","White color");
     sel->add_value(black,"black");
     sel->add_value(red,"red");
     sel->add_value(green,"green");
 }

 Record main("MainRecord", "The main record of flow.");
 main.declare_key("array_of_records", Array(output_record), "Array of output streams.");
 main.declare_key("record_record", record_record, "no comment on record_record");
 main.declare_key("color", sel, "My favourite color.");
 main.declare_key("color1", sel, "My second favourite color.");
 main.declare_key("array_record", array_record, "no commment on array_record");


 main.documentation(cout, true);
}
