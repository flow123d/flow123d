/**
 * type_output_test.cpp
 */

#include <flow_gtest.hh>

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

    Selection sel("Colors", "Selection of colors");
    {
        sel.add_value(blue, "blue");
        sel.add_value(white,"white","White color");
        sel.add_value(black,"black");
        sel.add_value(red,"red");
        sel.add_value(green,"green");
        sel.close();
    }

    Record main("MainRecord", "The main record of flow.");
    main.declare_key("array_of_records", Array(output_record), Default::obligatory(), "Array of output streams.");
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

    cout << "## " << "OutputJSONMachine printout" << endl;

    cout << OutputJSONMachine(&main) << endl;
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

    cout << "## " << "OutputText printout" << endl;
    OutputText output_text( &b_rec, 0);
    output_text.print(cout);

    OutputText output_text2( &c_rec, 0);
    output_text2.print(cout);

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(&a_rec) << endl;
}


/**
 * Child classes of Input::Type::AbstractRecord and AdHocAbstractRecord
 * Contains public method for adding descendants
 */
class AbstractRecordTest : public Input::Type::AbstractRecord {
public:
	AbstractRecordTest(const string & type_name_in, const string & description)
	: Input::Type::AbstractRecord(type_name_in, description)
	{}

	void declare_descendant(const Record &subrec) {
		add_descendant(subrec);
	}
};

TEST(OutputTypeAbstractRecord, ad_hoc_abstract_record_test) {
    using namespace Input::Type;

	Selection sel_problem("Problem_TYPE_selection");
	{
		sel_problem.add_value(0, "B_Record");
		sel_problem.add_value(1, "C_Record");
		sel_problem.add_value(2, "D_Record");
		sel_problem.add_value(3, "E_Record");
		sel_problem.close();
	}

	Record b_rec("B_Record", "Test record.");
	b_rec.declare_key("TYPE", sel_problem, Default("B_Record"), "Type of problem");
    b_rec.declare_key("b_val", Integer(), Default("10"), "");
    b_rec.declare_key("description", String(), Default::obligatory(), "");
    b_rec.close();

    Record c_rec("C_Record", "Test record.");
	c_rec.declare_key("TYPE", sel_problem, Default("C_Record"), "Type of problem");
    c_rec.declare_key("c_val", Double(), Default("0.5"), "");
    c_rec.declare_key("mesh", String(), Default("input.msh"), "Comp. mesh.");
    c_rec.close();

    Record d_rec("D_Record", "Test record.");
	d_rec.declare_key("TYPE", sel_problem, Default("D_Record"), "Type of problem");
    d_rec.declare_key("d_val", Integer(), Default("1"), "");
    d_rec.declare_key("pause", Bool(), Default("false"), "");
    d_rec.close();

    Record e_rec("E_Record", "Test record.");
    e_rec.declare_key("TYPE", sel_problem, Default("E_Record"), "Type of problem");
    e_rec.declare_key("e_val", String(), Default("Some value"), "");
    e_rec.declare_key("pause", Bool(), Default("false"), "");
    e_rec.close();

    // ancestor abstract record
    AbstractRecordTest a_rec_test("EqBase", "Base of equation records.");
    a_rec_test.close();
    a_rec_test.declare_descendant(b_rec);
    a_rec_test.declare_descendant(c_rec);
    AbstractRecord a_rec(a_rec_test);

    // adhoc abstract record - descendant of a_rec
    AdHocAbstractRecord adhoc_rec(a_rec);
    adhoc_rec.add_child(d_rec);
    adhoc_rec.add_child(e_rec);
    adhoc_rec.finish();

    Record root_rec("Root", "Root record.");
    root_rec.declare_key("problem", adhoc_rec, Default::obligatory(), "Base problem");
    root_rec.declare_key("pause", Bool(), Default("false"), "");
    root_rec.close();

    cout << "## " << "AdHocAbstractRecord OutputText printout";
    OutputText output_text( &root_rec, 0);
    output_text.print(cout);

    cout << endl << "## " << "OutputJSONTemplate printout";
    OutputJSONTemplate output_json_template( &root_rec, 0);
    output_json_template.print(cout);

    cout << endl << "## " << "OutputLatex printout";
    OutputLatex output_latex( &root_rec, 0);
    output_latex.print(cout);

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    OutputJSONMachine output_json_machine( &root_rec, 0);
    output_json_machine.print(cout);
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

#include <boost/regex.hpp>
TEST(OutputTypeRegEx, regex_filter_test) {
    //static const boost::regex e("(\\d{4}[- ]){3}\\d{4}");
    //regex_match("", e);

    using namespace Input::Type;

    Record a_rec("FieldConstant:Field:R3 -> Real", "");
    a_rec.close();

    // test derived type
    Record b_rec("FieldConstant:Field:R3 -> Real[3,3]", "");
    b_rec.close();

    Record c_rec("FieldConstant:Field:R3 -> Enum[3]", "");
    c_rec.close();


    Record main("MainRecord", "The main record of flow.");
    main.declare_key("a", a_rec, "first record of flow");
    main.declare_key("b", b_rec, "first record of flow");
    main.declare_key("c", c_rec, "second record of flow");
    main.close();

    OutputText output_text( &main, 0);
    // finds expressions in format '[N,N]' or '[N]' where N is dimension (0-3)
    //output_text.set_filter("(\\[[0-3]\\,?[0-3]?\\])");
    output_text.set_filter(":Field:.*");
    output_text.print(cout);
}
