/**
 * address_test.cpp
 */

#include <flow_gtest.hh>

#include "input/type_base.hh"
#include "input/accessors.hh"
#include "input/type_output.hh"
#include "input/reader_to_storage.hh"
#include "input/storage.hh"

namespace IT = Input::Type;

const string read_input_json = R"JSON(
{
	problem = {
		TYPE = "SequentialCoupling",
		regions = [
           { name = "region_1", id = 5, other={} },
           { name = "region_2", id = 10, other={} }
        ]
    }
}
)JSON";

DECLARE_EXCEPTION(ExcTest, << Input::EI_Address::val);

TEST(InputAddress, address_output_test) {
	IT::Selection sel_problem = IT::Selection("Problem_TYPE_selection")
		.add_value(0, "SequentialCoupling", "sequential coupling problem")
    	.add_value(1, "Other")
		.close();

	IT::Record other_record = IT::Record("Other", "Record with data for other problem")
			.declare_key("TYPE", IT::String(), IT::Default("Other"),	"Type of problem")
		    .close();
	{
		other_record.finish();
	}


	IT::Record region_input_type = IT::Record("Region", "Definition of region of elements.")
		.declare_key("name", IT::String(), IT::Default::obligatory(),
		                "Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
		                "The ID of the region to which you assign label.")
		.declare_key("other", other_record, IT::Default::obligatory(),
		                "The ID of the region to which you assign label.")
		.close();

	IT::Record sequential_coupling_record = IT::Record("SequentialCoupling", "Record with data for a general sequential coupling")
		.declare_key("TYPE", IT::String(), IT::Default("SequentialCoupling"),
				"Type of problem")
		.declare_key("regions", IT::Array(region_input_type),IT::Default::obligatory(),
	             "Definition of region.")
		.close();


    IT::Abstract problem = IT::Abstract("Problem","Base problem.")
    	.close();
    {
    	problem.add_child(sequential_coupling_record);
    	problem.add_child(other_record);
    }

    IT::Record root_record = IT::Record("RootRecord", "Root record of JSON input")
		.declare_key("problem", problem, IT::Default::obligatory(),
				"Simulation problem to be solved.")
		.close();

	Input::ReaderToStorage json_reader( read_input_json,  root_record, Input::FileFormat::format_JSON);
	Input::Record i_rec = json_reader.get_root_interface<Input::Record>();

	EXPECT_EQ("/", i_rec.address_string() );
	Input::AbstractRecord problem_rec = i_rec.val<Input::AbstractRecord>("problem");
	EXPECT_EQ("/problem", problem_rec.address_string() );
	Input::Array regions = Input::Record(problem_rec).val<Input::Array>("regions");
	EXPECT_EQ("/problem/regions", regions.address_string() );
	Input::Record reg = *(regions.begin<Input::Record>() );
	EXPECT_EQ("/problem/regions/0", reg.address_string() );
	EXPECT_EQ("/problem/regions/0/other", reg.val<Input::Record>("other").address_string() );

	EXPECT_THROW_WHAT({ THROW( ExcTest() << reg.ei_address());}, ExcTest, "Program Error: /problem/regions/0");

}

