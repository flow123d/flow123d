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

TEST(InputTypeArray, construction) {
using namespace Input::Type;
    Array arr_int(Integer());
    Array arr_arr_dbl( Array( Double() ));

    Record rec_1("record_type_1","desc");
    boost::shared_ptr<Record> rec_2 = boost::make_shared<Record>("record_type_2", "desc");


    ASSERT_DEATH( {Array arr_rec_ref( Record("subrec_type", "desc") ); },
                 "Complex type .* shared_ptr."
                );

    Array arr_rec_shared_ptr( rec_2 );
}

TEST(InputTypeSelection, declare_selection) {
using namespace Input::Type;


}

TEST(InputTypeRecord, declare_key_scalars) {
using namespace Input::Type;

   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FileName
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   static Record output_record("OutputRecord",
           "Information about one file for field data.");
   {
   // FileName, with default, reference, default read_time
   output_record.declare_key("file", FileName( output_file ), DefaultValue(DefaultValue::read_time),
           "File for output stream.");
   // Integer, with default, shared_ptr, dafault on declaration
   boost::shared_ptr<Integer> digits_type=boost::make_shared<Integer>((int)0, (int)8);
   output_record.declare_key("digits",digits_type, DefaultValue("8"),
           "Number of digits used for output double values into text output files.");
   // Bool, without default, reference
   output_record.declare_key("compression", Bool(),
           "Whether to use compression of otput file.");
   // Double, without dflt, shared_ptr
   boost::shared_ptr<Double> time_type=boost::make_shared<Double>(0.0);
   output_record.declare_key("start_time", time_type,
           "Simulation time of first output.");
   // String, with default == none, reference
   output_record.declare_key("data_description", String(),DefaultValue(),
           "");

   ASSERT_DEATH( { output_record.declare_key("data_description", String(),"");},
                "Redeclaration of key:");

   ASSERT_DEATH( { DefaultValue d(DefaultValue::declaration);},
                "Can not construct DefaultValue with type 'declaration' without providing the default value.");

   // should fail at compile time
   // output_record.declare_key("xx", 10, "desc");

   //boost::shared_ptr<int> int_ptr = boost::make_shared<int>(8);
   //output_record.declare_key("xx", int_ptr, "desc");


   } // delete local variables

   output_record.documentation(cout, true);
}

TEST(InputTypeRecord, declare_key_arrays) {
using namespace Input::Type;

   static Record array_record("RecordOfArrays",
            "Long description of record.\n"
            "Description could have more lines"
            );
   {
    // array type passed through shared_ptr
    boost::shared_ptr<Array> array_of_int = boost::make_shared<Array>(Integer(0), 5, 100 );
    array_record.declare_key("array_of_5_ints", array_of_int,
            "Some bizare array.");
    // array type passed by reference
    array_record.declare_key("array_of_str", Array( String() ), DefaultValue(),
            "Desc. of array");
    array_record.declare_key("array_of_str_1", Array( String() ), DefaultValue(),
                "Desc. of array");
   }

    array_record.documentation(cout, true);

    ASSERT_DEATH( {array_record.declare_key("some_key", Array( String() ), DefaultValue("10"), ""); },
                  "Can not provide default value for non scalar type."
                 );

}

TEST(InputTypeRecord, declare_key_record) {
using namespace Input::Type;

    Record record_record("RecordOfRecords",
            "Long description of record.\n"
            "Description could have more lines"
            );


    {
        ASSERT_DEATH( {record_record.declare_key("sub_rec_1", Record("subrec_type", "desc") , "desc"); },
                     "Complex type .* shared_ptr."
                    );


        boost::shared_ptr<Record> other_record=boost::make_shared<Record>("OtherRecord","desc");
        record_record.declare_key("sub_rec_1", other_record, "key desc");
    }

    record_record.documentation(cout, true);

}

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

