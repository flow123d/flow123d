/**
 * reference_test.cpp
 */

#include <flow_gtest.hh>

#include "input/type_base.hh"
#include "input/accessors.hh"
#include "input/type_output.hh"
#include "input/json_to_storage.hh"

namespace IT = Input::Type;

const string valid_record_json = R"JSON(
{
	data = {
		name = "Some name",
        ids = [0, 5, 12 ],
        description = "Some text",
        value = {
            some_boolean = false,
            val_description = {REF="/data/description"}
        }
	},
	pause_after_run = false
}
)JSON";

const string cyclic_record_json = R"JSON(
{
	data = {
		name = "Some name",
        ids = [0, 5, 12 ],
        description = "Some text",
        value = {REF="/data/value"}
	},
	pause_after_run = false
}
)JSON";

const string cyclic_array_json = R"JSON(
{
	data = {
		name = "Some name",
        ids = [0, 5, {REF="/data/ids/2"} ],
        description = "Some text",
        value = {
            some_boolean = false,
            val_description = "Short description"
        }
	},
	pause_after_run = false
}
)JSON";


IT::Record get_type_record() {
    IT::Record value_rec("value", "Record with boolean value");
    {
    	value_rec.declare_key("some_boolean", IT::Bool(), IT::Default("false"),
             "Some boolean value.");
    	value_rec.declare_key("val_description", IT::String(), IT::Default::optional(),
             "Description of value.");
    	value_rec.close();
    }

	IT::Record data_rec("Data", "Record with data");
	{
		data_rec.declare_key("name", IT::String(), IT::Default::obligatory(),
	             "Name of data part.");
		data_rec.declare_key("ids", IT::Array( IT::Integer(0) ), IT::Default::optional(),
		                "The IDs of the value.");
		data_rec.declare_key("description", IT::String(), IT::Default::optional(),
	             "Short description.");
		data_rec.declare_key("value", value_rec,
	             "Value record.");
		data_rec.close();
	}

	IT::Record root_record("Problem", "Record of problem");
	{
		root_record.declare_key("data", data_rec, IT::Default::obligatory(),
	             "Definition of data.");
    	root_record.declare_key("pause_after_run", IT::Bool(), IT::Default("false"),
             "If true, the program will wait for key press before it terminates.");
		root_record.close();
	}

	return root_record;
}



TEST(JSONReference, valid_reference_rec_test) {
	using namespace Input;

    JSONToStorage json_reader;
	stringstream ss(valid_record_json.c_str());
	json_reader.read_stream( ss, get_type_record());
	// Input::Record i_rec = json_reader.get_root_interface<Input::Record>();
}

TEST(JSONReference, cyclic_reference_rec_test) {
	using namespace Input;

    JSONToStorage json_reader;
	stringstream ss(cyclic_record_json.c_str());
    EXPECT_THROW_WHAT( {json_reader.read_stream( ss, get_type_record());}, JSONPath::ExcReferenceNotFound,
			"cannot follow reference");
}

TEST(JSONReference, cyclic_reference_arr_test) {
	using namespace Input;

    JSONToStorage json_reader;
	stringstream ss(cyclic_array_json.c_str());
    EXPECT_THROW_WHAT( {json_reader.read_stream( ss, get_type_record());}, JSONPath::ExcReferenceNotFound,
			"cannot follow reference");
}
