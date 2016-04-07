/**
 * type_output_test.cpp
 */

#include <flow_gtest.hh>

#include "input/type_base.hh"
#include "input/type_output.hh"
#include "input/type_generic.hh"
#include "input/type_tuple.hh"

// Set version of program and other base data. Simulate data got from rev_num.h.
Input::Type::RevNumData get_rev_num_data() {
	static Input::Type::RevNumData rev_num_data;
	rev_num_data.version = "X.X.X";
	rev_num_data.revision = "unit.test";
	rev_num_data.branch = "unit.test";
	rev_num_data.url = "unit.test";

	return rev_num_data;
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

    Tuple tuple("TimeStepping",
            "Basic time information.");
    {
    	tuple.declare_key("start_time", Double(0), Default::obligatory(), "Start time.");
    	tuple.declare_key("end_time", Double(0), Default::obligatory(), "End time.");
    	tuple.declare_key("time_step", Double(0), Default("0.1"), "Time step.");
        tuple.close();
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
    main.declare_key("time", tuple, "Time stepping management");
    main.close();

    cout << "## " << "OutputText printout" << endl;
    cout << OutputText(&main) << endl << endl;

    cout << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(get_rev_num_data()) << endl;
}

TEST(OutputTypeAbstract, abstract_record_test) {
    using namespace Input::Type;

	Record copy_rec = Record("Copy","")
       	.declare_key("mesh", String(), Default("\"input.msh\""), "Comp. mesh.")
       	.declare_key("a_val", String(), Default::obligatory(), "")
		.close();

	Abstract a_rec = Abstract("EqBase","Base of equation records.");
	Abstract &a_ref = a_rec.allow_auto_conversion("EqDarcy").close();
    EXPECT_EQ( a_rec, a_ref);

    // test derived type
    Record b_rec = Record("EqDarcy","test derived type and reducible key")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
    	.declare_key("b_val", Integer(), Default("10"), "")
    	.allow_auto_conversion("a_val")
		.close();

    Record c_rec = Record("EqTransp","test derived type")
    	.derive_from(a_rec)
		.copy_keys(copy_rec)
    	.declare_key("c_val", Integer(), "")
    	.declare_key("a_val", Double(),"")
		.close();

    cout << "## " << "OutputText printout" << endl;
    cout << OutputText( &b_rec);
    cout << OutputText( &c_rec);

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(get_rev_num_data()) << endl;
}


TEST(OutputTypeAbstract, ad_hoc_abstract_test) {
    using namespace Input::Type;

    // ancestor abstract record
	Abstract a_rec = Abstract("EquationBase", "Base of equation records.")
    		.close();

    Record b_rec = Record("B_Record", "Test record.")
    	.derive_from(a_rec)
    	.declare_key("b_val", Integer(), Default("10"), "")
    	.declare_key("description", String(), Default::obligatory(), "")
    	.close();

    Record c_rec = Record("C_Record", "Test record.")
		.derive_from(a_rec)
    	.declare_key("c_val", Double(), Default("0.5"), "")
    	.declare_key("mesh", String(), Default("\"input.msh\""), "Comp. mesh.")
    	.close();

    Record d_rec = Record("D_Record", "Test record.")
		.declare_key("TYPE", String(), Default("\"D_Record\""), "Type of problem")
    	.declare_key("d_val", Integer(), Default("1"), "")
    	.declare_key("pause", Bool(), Default("false"), "")
    	.close();

    Record e_rec = Record("E_Record", "Test record.")
    	.declare_key("TYPE", String(), Default("\"E_Record\""), "Type of problem")
    	.declare_key("e_val", String(), Default("\"Some value\""), "")
    	.declare_key("pause", Bool(), Default("false"), "")
    	.close();

    // adhoc abstract - descendant of a_rec
    AdHocAbstract adhoc_rec = AdHocAbstract(a_rec).close();
    adhoc_rec.add_child(d_rec);
    adhoc_rec.add_child(e_rec);
    adhoc_rec.finish();

    Record root_rec = Record("Root", "Root record.")
    	.declare_key("problem", adhoc_rec, Default::obligatory(), "Base problem")
    	.declare_key("pause", Bool(), Default("false"), "")
    	.close();

    cout << "## " << "OutputText printout";
    cout << OutputText( &root_rec) << endl;

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(get_rev_num_data()) << endl;
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
    cout << OutputText(&array_record) << endl;

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(get_rev_num_data()) << endl;
}


TEST(OutputTypeParameter, parameter_test) {
    using namespace Input::Type;

    std::vector<TypeBase::ParameterPair> param_vec;
    param_vec.push_back( std::make_pair("param", boost::make_shared<Integer>()) );

    static Record param_record = Record("WithParameter", "Record with parameter.")
			.declare_key("param", Parameter("param"), "desc.")
			.declare_key("start_time", Double(), "desc.")
			.declare_key("name", String(), "desc.")
			.close();

	static Instance inst = Instance(param_record, param_vec)
								.close();

    cout << endl << "## " << "OutputJSONMachine printout" << endl;
    cout << OutputJSONMachine(get_rev_num_data()) << endl;
}
