/*
 * type_record_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */


#include <gtest/gtest.h>
#include "gtest_throw_what.hh"

#include <input/type_record.hh>



/**
 * Test Record class.
 */
TEST(InputTypeRecord, declare_key_scalars) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FilePath
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   static Record rec("SomeRecord", "desc.");

   rec.declare_key("file", FileName::output(), Default::optional(), "desc.");

   Integer digits_type(0, 8);
   rec.declare_key("digits",digits_type, Default("8"), "desc.");

   rec.declare_key("compression", Bool(),"desc.");

   Double time_type(0.0);
   rec.declare_key("start_time", time_type,"desc.");

   rec.declare_key("data_description", String(), Default::optional(),"");


   // errors during declaration
   Record rec_empty;
   EXPECT_DEATH( {rec_empty.declare_key("xx", Integer(), "");}, "Empty Record handle.");

   Record rec_fin("xx","");
   rec_fin.finish();
   EXPECT_DEATH( {rec_fin.declare_key("xx", String(),"");}, "in finished Record type:");

   EXPECT_DEATH( {rec.declare_key("ar",Array(Integer()), Default("[0, 1]"), "");} , "Default value for non scalar type in declaration of key:");

   Record rec_unfin("yy","");
   EXPECT_DEATH({rec.declare_key("yy", rec_unfin, ""); }, "Unfinished type of declaring key:");

   EXPECT_DEATH( { rec.declare_key("data_description", String(),"");},
                "Re-declaration of the key:");

   EXPECT_THROW_WHAT( { rec.declare_key("wrong_double", Double(), Default("1.23 4"),"");}, ExcWrongDefault,
           "Default value .* do not match type: 'Double';");

   enum Colors {
       white, black, red
   };

   Selection sel("Color selection");
   sel.add_value(black, "black");
   sel.add_value(red, "red");
   sel.finish();

   rec.declare_key("plot_color", sel, "Color to plot the fields in file.");

   // test correct finishing.
   EXPECT_DEATH( {rec.size();}, "Asking for information of unfinished Record type");

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
    Array array_of_int(Integer(0), 5, 100 );
    array_record.declare_key("array_of_5_ints", array_of_int,"Some bizare array.");

    // array type passed by reference
    array_record.declare_key("array_of_str", Array( String() ),"Desc. of array");
    array_record.declare_key("array_of_str_1", Array( String() ), "Desc. of array");


    ASSERT_DEATH( {array_record.declare_key("some_key", Array( String() ), Default("10"), ""); },
                  "Default value for non scalar type in declaration of key:"
                 );
    array_record.finish();
}

TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";


    Record record_record("RecordOfRecords", "");

    // Test that Record has to be passed as shared_ptr
    //ASSERT_DEATH( {record_record->declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
    //              "Complex type .* shared_ptr."
    //              );

    Record other_record("OtherRecord","desc");
    other_record.finish();

    record_record.declare_key("sub_rec_1", other_record, "key desc");

    // recursion  -  forbidden
    //record_record->declare_key("sub_rec_2", record_record, "desc.");

    record_record.finish();
}

TEST(InputTypeRecord, iterating) {
    using namespace Input::Type;

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    {
        output_record.declare_key("file", FileName::output(), Default::optional(),
                "File for output stream.");

        Integer digits_type((int)0, (int)8);
        output_record.declare_key("digits",digits_type, Default("8"),
                "Number of digits used for output double values into text output files.");
        output_record.declare_key("compression", Bool(),
                "Whether to use compression of output file.");

        Double time_type(0.0);
        output_record.declare_key("start_time", time_type,
                "Simulation time of first output.");
        output_record.declare_key("data_description", String(),
                "");
        output_record.finish();
    } // delete local variables

    // methods begin() and end(), basic work with iterators
    Record::KeyIter it = output_record.begin();
    EXPECT_EQ( 0, it->key_index);
    EXPECT_EQ("file", it->key_);
    EXPECT_EQ("File for output stream.", it->description_);
    EXPECT_EQ(typeid(FileName) , typeid(*(it->type_)));
    EXPECT_EQ(FilePath::output_file, static_cast<const FileName *>( &(*it->type_) )->get_file_type() );
    it+=4;
    EXPECT_EQ( "data_description", it->key_);
    EXPECT_EQ( output_record.end(), it+1 );
    // method size()
    EXPECT_EQ(5, output_record.size());
    //method key_index
    EXPECT_EQ(0,output_record.key_index("file"));
    EXPECT_EQ(2,output_record.key_index("compression"));
    EXPECT_THROW({output_record.key_index("x_file");}, Record::ExcRecordKeyNotFound );
    // method key_iterator
    EXPECT_EQ(2,output_record.key_iterator("compression")->key_index);


}

TEST(InputTypeRecord, check_key_validity) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord",
            "Information about one file for field data.");
    EXPECT_DEATH( {output_record.declare_key("a b",Bool(),"desc."); },
            "Invalid key identifier"
            );
    EXPECT_DEATH( {output_record.declare_key("AB",Bool(),"desc."); },
            "Invalid key identifier"
            );
    EXPECT_DEATH( {output_record.declare_key("%$a",Bool(),"desc."); },
            "Invalid key identifier"
            );

}

TEST(InputTypeRecord, RecordCopy) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    Record output_record("OutputRecord", "");
    output_record.declare_key("file", FileName::output(), "");


    Record copy_rec = output_record;

    Integer digits_type((int)0, (int)8);
    copy_rec.declare_key("digits",digits_type, "");
    copy_rec.finish();

    EXPECT_EQ( true, output_record.is_finished());
    EXPECT_EQ( true, output_record.has_key("digits") );
}

/**
 * Test Abstract Record.
 */

TEST(InputTypeAbstractRecord, inheritance) {
using namespace Input::Type;
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    AbstractRecord a_rec("EqBase","Base of equation records.");
    a_rec.declare_key("mesh", String(), Default::obligatory(), "Comp. mesh.");
    a_rec.declare_key("a_val", String(), Default::obligatory(), "");
    a_rec.finish();

    EXPECT_EQ(0, a_rec.key_index("TYPE"));
    EXPECT_EQ(Selection("EqBase_selection"), *(a_rec.key_iterator("TYPE")->type_ ));


    Record b_rec("EqDarcy","");
    b_rec.derive_from(a_rec);
    b_rec.declare_key("b_val", Integer(), "");
    b_rec.finish();
    EXPECT_EQ(0, b_rec.key_index("TYPE"));
    EXPECT_EQ(a_rec.key_iterator("TYPE")->type_, b_rec.key_iterator("TYPE")->type_);

    Record c_rec("EqTransp","");
    c_rec.derive_from(a_rec);
    c_rec.declare_key("c_val", Integer(), "");
    c_rec.declare_key("a_val", Double(),"");
    c_rec.finish();
    EXPECT_EQ(0, c_rec.key_index("TYPE"));

    a_rec.no_more_descendants();

    // inherited keys
    EXPECT_TRUE( b_rec.has_key("mesh") );
    EXPECT_TRUE( c_rec.has_key("mesh") );
    // overwritten key
    EXPECT_EQ( Double(), *(c_rec.key_iterator("a_val")->type_));

    //get descendant
    EXPECT_EQ( b_rec, a_rec.get_descendant("EqDarcy"));
    EXPECT_EQ( c_rec, a_rec.get_descendant("EqTransp"));

    EXPECT_THROW_WHAT( {c_rec.derive_from(a_rec);} , Record::ExcDeriveNonEmpty,
            "Can not derive from Record .* into non-empty Record:" );

    AbstractRecord x_rec("ar","");
    Record y_rec("y_rec","");
    EXPECT_DEATH( {y_rec.derive_from(x_rec);} , "Can not add descendant to unfinished AbstractType." );
}


