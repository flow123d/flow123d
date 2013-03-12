/**
 * type_output_test.cpp
 */

#include <gtest/gtest.h>

#include "input/type_base.hh"
#include "input/type_output.hh"

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

TEST(OutputTypeTypeBase, record_output_test) {
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

 cout << "## " << "OutputText printout" << endl;

 OutputText output_text( &main, 0);
 output_text.print(cout);

 cout << endl;
 cout << "## " << "OutputJSONTemplate printout" << endl;

 OutputJSONTemplate output_json( &main, 0);
 output_json.print(cout);

 cout << endl;
 cout << "## " << "OutputLatex printout" << endl;

 cout << OutputLatex(&main) << endl;


}

TEST(OutputTypeAbstractRecord, abstract_record_test) {
	using namespace Input::Type;

    AbstractRecord a_rec("EqBase","Base of equation records.");
    a_rec.declare_key("mesh", String(), Default("input.msh"), "Comp. mesh.");
    a_rec.declare_key("a_val", String(), Default::obligatory(), "");
    AbstractRecord &a_ref = a_rec.allow_auto_conversion("EqDarcy");
    EXPECT_EQ( a_rec, a_ref);
    a_rec.close();

    // test derived type
    Record b_rec("EqDarcy", "test derived type and reducible key");
    b_rec.derive_from(a_rec);
    b_rec.declare_key("b_val", Integer(), Default("10"), "");
    b_rec.allow_auto_conversion("a_val");

    Record c_rec("EqTransp","test derived type");
    c_rec.derive_from(a_rec);
    c_rec.declare_key("c_val", Integer(), "");
    c_rec.declare_key("a_val", Double(),"");

    c_rec.close();
    b_rec.close();

    OutputText output_text( &b_rec, 0);
    output_text.print(cout);

    OutputText output_text2( &c_rec, 0);
    output_text2.print(cout);

}

TEST(OutputTypeArray, array_of_array_test) {
	using namespace Input::Type;

	Record array_record("RecordOfArrays",
	         "Record contains array of array"
	         );
	{
		Array array_of_int(Integer(0), 0, 100 );
		Array array_of_array(array_of_int, 0, 100 );

		array_record.declare_key("array_of_array", array_of_array,
				 "Matrix of integer.");
		array_record.declare_key("data_description", String(), Default::optional(),
	            "");
		array_record.close();
	}

	cout << "## " << "OutputText printout" << endl;
    OutputText output_text( &array_record, 0);
    output_text.print(cout);

    cout << endl;
    cout << "## " << "OutputJSONTemplate printout" << endl;
    OutputJSONTemplate output_json( &array_record, 0);
    output_json.print(cout);

}
