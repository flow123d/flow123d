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
		description = "Basic solved problem.",
		mesh_file = "./input.msh",
		regions = {
            name = "region_1",
            id = 5,
            element_list = [ 0, 3, 5, 12 ]
        } 
	},
	pause_after_run = false
}
)JSON";

TEST(InputAddress, address_output_test) {
	IT::Selection sel_problem("Problem_TYPE_selection");
	{
		sel_problem.add_value(0, "SequentialCoupling", "sequential coupling problem");
		sel_problem.add_value(1, "Other");
		sel_problem.close();
	}

	IT::Record region_input_type("Region", "Definition of region of elements.");
	{
		region_input_type.declare_key("name", IT::String(), IT::Default::obligatory(),
		                "Label (name) of the region. Has to be unique in one mesh.\n");
		region_input_type.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
		                "The ID of the region to which you assign label.");
		region_input_type.declare_key("element_list", IT::Array( IT::Integer(0) ), IT::Default::optional(),
		                "Specification of the region by the list of elements. This is not recomended");
		region_input_type.close();
	}

	IT::Record sequential_coupling_record("SequentialCoupling", "Record with data for a general sequential coupling");
	{
		sequential_coupling_record.declare_key("TYPE", sel_problem, IT::Default("SequentialCoupling"),
				"Type of problem");
		sequential_coupling_record.declare_key("description", IT::String(), IT::Default::optional(),
	             "Short description of the solved problem.");
		sequential_coupling_record.declare_key("mesh_file", IT::FileName::input(), IT::Default::obligatory(),
	             "Input file with mesh description.");
		sequential_coupling_record.declare_key("regions", region_input_type,
	             "Definition of region.");
		sequential_coupling_record.close();
	}

	IT::Record other_record("Other", "Record with data for other problem");
	{
		other_record.declare_key("TYPE", sel_problem, IT::Default("Other"),
				"Type of problem");
		other_record.declare_key("description", IT::String(), IT::Default::optional(),
	             "Short description of the solved problem.");
		other_record.declare_key("input_file", IT::FileName::input(), IT::Default::obligatory(),
	             "Input file of solved problem.");
		other_record.declare_key("steps", IT::Integer(0),
	             "Count of steps.");
		other_record.close();
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
    	root_record.declare_key("pause_after_run", IT::Bool(), IT::Default("false"),
             "If true, the program will wait for key press before it terminates.");
    	root_record.close();
    }

    Input::JSONToStorage json_reader;
	stringstream ss(read_input_json.c_str());
	json_reader.read_stream( ss,  root_record);
	Input::Record i_rec = json_reader.get_root_interface<Input::Record>();

	Input::Address a_problem( *(i_rec.get_address().down(0)) );
	Input::Address a_description( *(a_problem.down(1)) );
	Input::Address a_regions( *(a_problem.down(3)) );
	Input::Address a_element_list( *(a_regions.down(2)) );
	Input::Address a_element_1( *(a_element_list.down(1)) );

	EXPECT_EQ("/", i_rec.get_address().make_full_address());
	EXPECT_EQ("/problem", a_problem.make_full_address());
	EXPECT_EQ("/problem/description", a_description.make_full_address());
	EXPECT_EQ("/problem/regions", a_regions.make_full_address());
	EXPECT_EQ("/problem/regions/element_list", a_element_list.make_full_address());
	EXPECT_EQ("/problem/regions/element_list/1", a_element_1.make_full_address());
}
