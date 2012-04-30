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
    rec_2->finish();

    //MAKE_Input_Type_Record(tmp)("subrec_type", "desc");
    //EXPECT_DEATH( {Array arr_rec_ref( *tmp ); },
    //             "Complex type .* shared_ptr."
    //            );

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
    green,
    yellow
};

TEST(InputTypeSelection, construction) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Selection<enum Colors> sel("Colors");
    sel.add_value(blue, "blue");
    sel.add_value(white,"white","White color");
    sel.add_value(black,"black");
    sel.add_value(red,"red");
    EXPECT_DEATH( {sel.add_value(green,"red");}, "already exists in Selection:");
    EXPECT_DEATH( {sel.add_value(blue,"blue1");}, "conflicts with value");
    sel.add_value(green,"green");
    sel.finish();

    // getter methods
    EXPECT_TRUE( sel.has_name("blue") );
    EXPECT_FALSE( sel.has_name("xblue") );
    EXPECT_TRUE( sel.has_value(blue) );
    EXPECT_TRUE( sel.has_value(black) );
    EXPECT_FALSE( sel.has_value(yellow) );

    EXPECT_EQ( 45, sel.name_to_value("black") );
    EXPECT_THROW( {sel.name_to_value("xblack");}, SelectionKeyNotFound );


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

   rec.declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::optional), "desc.");

   boost::shared_ptr<Integer> digits_type=boost::make_shared<Integer>((int)0, (int)8);
   rec.declare_key("digits",digits_type, DefaultValue("8"), "desc.");

   rec.declare_key("compression", Bool(),"desc.");

   boost::shared_ptr<Double> time_type=boost::make_shared<Double>(0.0);
   rec.declare_key("start_time", time_type,"desc.");

   rec.declare_key("data_description", String(),DefaultValue(),"");



   EXPECT_DEATH( { rec.declare_key("data_description", String(),"");},
                "Re-declaration of the key:");

   EXPECT_DEATH( { DefaultValue d(DefaultValue::declaration);},
                "Can not construct DefaultValue with type 'declaration' without providing the default value.");

   enum Colors {
       white, black, red
   };

   MAKE_Input_Type_Selection(Colors, sel)("Color selection");
   sel->add_value(black, "black");
   sel->add_value(red, "red");
   sel->finish();

   rec.declare_key("plot_color", sel, "Color to plot the fields in file.");

   // test correct finishing.
   EXPECT_DEATH( {rec.documentation(cout);}, "Can not provide documentation of unfinished Record type: ");
   rec.finish();
   EXPECT_DEATH( {rec.declare_key("xx", String(),"");}, "in finished Record type:");

   //EXPECT_DEATH( { rec.declare_key("x", *sel, "desc.");},
   //             "Complex type .* shared_ptr.");

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
                  "Default value for non scalar type in declaration of key:"
                 );
    array_record.finish();
}

TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


    boost::shared_ptr<Record> record_record=boost::make_shared<Record>("RecordOfRecords", "");

    // Test that Record has to be passed as shared_ptr
    //ASSERT_DEATH( {record_record->declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
    //              "Complex type .* shared_ptr."
    //              );

    boost::shared_ptr<Record> other_record=boost::make_shared<Record>("OtherRecord","desc");
    other_record->finish();

    record_record->declare_key("sub_rec_1", other_record, "key desc");

    // direct recursion (indirect is forbidden)
    record_record->declare_key("sub_rec_2", record_record, "desc.");

    record_record->finish();
}

TEST(InputTypeRecord, iterating) {
    using namespace Input::Type;

    boost::shared_ptr<Record> output_record=boost::make_shared<Record>("OutputRecord",
            "Information about one file for field data.");
    {
        output_record->declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::optional),
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
        output_record->finish();
    } // delete local variables

    // methods begin() and end(), basic work with iterators
    Record::KeyIter it = output_record->begin();
    EXPECT_EQ( 0, it->key_index);
    EXPECT_EQ("file", it->key_);
    EXPECT_EQ("File for output stream.", it->description_);
    EXPECT_EQ(typeid(FileName) , typeid(*(it->type_)));
    EXPECT_EQ(output_file, static_cast<const FileName *>( &(*it->type_) )->get_file_type() );
    it+=4;
    EXPECT_EQ( "data_description", it->key_);
    EXPECT_EQ( output_record->end(), it+1 );
    // method size()
    EXPECT_EQ(5, output_record->size());
    //method key_index
    EXPECT_EQ(0,output_record->key_index("file"));
    EXPECT_EQ(2,output_record->key_index("compression"));
    EXPECT_THROW({output_record->key_index("x_file");}, KeyNotFound );
    // method key_iterator
    EXPECT_EQ(2,output_record->key_iterator("compression")->key_index);


}

TEST(InputTypeRecord, check_key_validity) {
    using namespace Input::Type;

    boost::shared_ptr<Record> output_record=boost::make_shared<Record>("OutputRecord",
            "Information about one file for field data.");
    EXPECT_DEATH( {output_record->declare_key("a b",Bool(),"desc."); },
            "Invalid key identifier"
            );
    EXPECT_DEATH( {output_record->declare_key("AB",Bool(),"desc."); },
            "Invalid key identifier"
            );
    EXPECT_DEATH( {output_record->declare_key("%$a",Bool(),"desc."); },
            "Invalid key identifier"
            );

}
/**
 * Test documentation output.
 */
TEST(InputTypeDocumentation, whole_tree) {
using namespace Input::Type;

boost::shared_ptr<Record> output_record=boost::make_shared<Record>("OutputRecord",
        "Information about one file for field data.");
{
    output_record->declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::optional),
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
    output_record->finish();
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
 array_record->finish();
}


 boost::shared_ptr<Record> record_record= boost::make_shared<Record>("RecordOfRecords",
         "Long description of record.\n"
         "Description could have more lines"
         );


 {
     boost::shared_ptr<Record> other_record=boost::make_shared<Record>("OtherRecord","desc");
     other_record->finish();

     record_record->declare_key("sub_rec_1", other_record, "key desc");

     // recursion
     record_record->declare_key("sub_rec_2", record_record, "Recursive key.");

     record_record->finish();
 }

 MAKE_Input_Type_Selection(enum Colors, sel)("Colors");
 {
     sel->add_value(blue, "blue");
     sel->add_value(white,"white","White color");
     sel->add_value(black,"black");
     sel->add_value(red,"red");
     sel->add_value(green,"green");
     sel->finish();
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
