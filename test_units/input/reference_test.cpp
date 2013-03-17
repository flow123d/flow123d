/**
 * reference_test.cpp
 */

#include <gtest/gtest.h>

#include "input/type_base.hh"
#include "input/accessors.hh"
#include "input/type_output.hh"
#include "input/json_to_storage.hh"

namespace IT = Input::Type;

const string read_input_json = R"JSON(
{
	data = {
		name = "Some name",
        ids = [0, 5, 12 ],
        description = "Some text",
        value = {REF="/data"}
	},
	pause_after_run = false
}
)JSON";

TEST(InputStorage, reference_test) {
    IT::Record value_rec("value", "Test record, not in JSON");
    {
    	value_rec.declare_key("some_boolean", IT::Bool(), IT::Default("false"),
             "Some boolean value.");
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


    Input::JSONToStorage json_reader;
	stringstream ss(read_input_json.c_str());
	EXPECT_DEATH( { json_reader.read_stream( ss,  root_record); },
			"JSON input contains cyclic reference:");
	// Input::Record i_rec = json_reader.get_root_interface<Input::Record>();
}
