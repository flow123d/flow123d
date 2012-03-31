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
 * - testovat zda pri declare_key funguje staticky assert na typ vs. defaultni hodnota vs. zadani shared_ptr resp. reference
 */

TEST(input_type_test, read) {
using namespace Input::Type;

   // make auxiliary record and test declare_key for
   // - various Scalar types (excluding Selection): Integer, Bool, Double, String, FileName
   // - various decalre_key templates: with/without default, shared_ptr/ reference
   // - various default values
   static Record output_record("OutputRecord",
           "Information about one file for field data.");
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


   static Record main_record("MainRecord",
            "Main application record.\n"
            "Description could have more lines"
            );

    Integer int_type(0,3);
    main_record.declare_key("size", int_type,  DefaultValue("2"),
            "Size of element array.");
    main_record.declare_key("tolerance", Double(0,0.001),  DefaultValue("2.3E-6"),
            "Tolerance of solver.");
    main_record.declare_key("substances", Array( String() ),
            "Substances names." );


    main_record.documentation(cout, true);

}
