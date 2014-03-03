/**
 * address_test.cpp
 */

#include <flow_gtest.hh>

#include "input/type_base.hh"
#include "input/accessors.hh"
#include "input/type_output.hh"
#include "input/json_to_storage.hh"
#include "input/storage.hh"

namespace IT = Input::Type;

/**
 * Child class of Input::Type::AbstractRecord
 * Contains public method for adding descendants
 */
class AbstractRecordTest : public IT::AbstractRecord {
public:
	AbstractRecordTest(const string & type_name_in, const string & description) : IT::AbstractRecord(type_name_in, description)
	{}

	void declare_descendant(const Record &subrec) {
		add_descendant(subrec);
	}
};

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
	IT::Selection sel_problem("Problem_TYPE_selection");
	sel_problem.add_value(0, "SequentialCoupling", "sequential coupling problem");
    sel_problem.add_value(1, "Other");
	sel_problem.close();

	IT::Record other_record = IT::Record("Other", "Record with data for other problem")
			.declare_key("TYPE", sel_problem, IT::Default("Other"),	"Type of problem")
		    .close();


	IT::Record region_input_type("Region", "Definition of region of elements.");
	{
		region_input_type.declare_key("name", IT::String(), IT::Default::obligatory(),
		                "Label (name) of the region. Has to be unique in one mesh.\n");
		region_input_type.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
		                "The ID of the region to which you assign label.");
		region_input_type.declare_key("other", other_record, IT::Default::obligatory(),
		                "The ID of the region to which you assign label.");
		region_input_type.close();
	}

	IT::Record sequential_coupling_record("SequentialCoupling", "Record with data for a general sequential coupling");
	{
		sequential_coupling_record.declare_key("TYPE", sel_problem, IT::Default("SequentialCoupling"),
				"Type of problem");
		sequential_coupling_record.declare_key("regions", IT::Array(region_input_type),IT::Default::obligatory(),
	             "Definition of region.");
		sequential_coupling_record.close();
	}


	AbstractRecordTest problemTest("Problem","Base problem.");
    {
		problemTest.close();
		problemTest.declare_descendant(sequential_coupling_record);
		problemTest.declare_descendant(other_record);
    }

    IT::AbstractRecord problem(problemTest);

    IT::Record root_record("RootRecord", "Root record of JSON input");
    {
    	root_record.declare_key("problem", problem, IT::Default::obligatory(),
             "Simulation problem to be solved.");
    }

    Input::JSONToStorage json_reader;
	stringstream ss(read_input_json.c_str());
	json_reader.read_stream( ss,  root_record);
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

