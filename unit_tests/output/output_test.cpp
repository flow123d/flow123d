/*
 * vtk_output_test.cpp
 *
 *      Author: Jiri Hnidek <jiri.hnidek@tul.cz>
 */

#define TEST_USE_MPI
#include <gtest_mpi.hh>

#include "io/output.h"

#include "input/json_to_storage.hh"
#include "input/accessors.hh"

#include "system/sys_profiler.hh"

#include "fields/field_base.hh"
#include "fields/field_constant.hh"

// Test input for output stream
const string output_stream = R"JSON(
{
  file = "./test1.pvd", 
  format = {
    TYPE = "vtk", 
    variant = "ascii"
  }, 
  name = "flow_output_stream"
}
)JSON";


// Test input for output data
const string foo_output = R"JSON(
{
  pressure_p0 = "flow_output_stream",
  material_id = "flow_output_stream"
}
)JSON";

using namespace Input::Type;

class Foo {
public:
    Foo() {}
    ~Foo() {}
    static Record input_type;
};

Record Foo::input_type
    = Record("Foo", "Output setting for foo equations.")
    .declare_key("pressure_p0", String(),
            "Name of output stream for P0 approximation of pressure.")
    .declare_key("material_id", String(),
        "Name of output stream for material ID.");

class OutputTest : public testing::Test, public OutputTime {
protected:
    virtual void SetUp() {}
    virtual void TearDown() {
    	OutputTime::destroy_all();
    }
};

TEST( OutputTest, test_create_output_stream ) {
    Input::JSONToStorage reader_output_stream;

    /* Read input for output stream */
    std::stringstream in_output_stream(output_stream);
    reader_output_stream.read_stream(in_output_stream, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream.get_root_interface<Input::Record>());

    /* Make sure that there is one OutputTime instance */
    ASSERT_EQ(OutputTime::output_streams.size(), 1);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = (*output_iter);

    /* Pointer at OutputTime could not be NULL. Otherwise other assertation
     * would crash */
    ASSERT_TRUE(output_time != NULL);

    /* The name of instance has to be equal to configuration file:
     * "variable output_stream" */
    EXPECT_EQ(*(output_time->name), "flow_output_stream");

    /* The type of instance should be OutputVTK */
    EXPECT_EQ(output_time->file_format, OutputTime::VTK);

    /* The filename should be "./test1.pvd" */
    EXPECT_EQ(*(output_time->base_filename()), "./test1.pvd");
}


TEST( OutputTest, find_outputstream_by_name ) {
    /* There should be at least one output stream */
    ASSERT_EQ(OutputTime::output_streams.size(), 1);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream"));
}


TEST( OutputTest, test_register_elem_fields_data ) {
    Input::JSONToStorage reader_output;
    /* Read input for output */
    std::stringstream in_output(foo_output);
    reader_output.read_stream(in_output, Foo::input_type);

    Profiler::initialize();

    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Field<3, FieldValue<1>::Scalar> scalar_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    scalar_field.set_default( Input::Type::Default("2") );
    scalar_field.set_name("pressure_p0");
    scalar_field.set_units("L");
    scalar_field.set_mesh(&mesh);
    scalar_field.set_time(0.0);

    /* Register scalar (double) data */
    OutputTime::register_data<3, FieldValue<1>::Scalar>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::ELEM_DATA, &scalar_field);

    Field<3, FieldValue<1>::Integer> integer_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    integer_field.set_default( Input::Type::Default("10") );
    integer_field.set_name("material_id");
    integer_field.set_units("");
    integer_field.set_mesh(&mesh);
    integer_field.set_time(0.0);

    /* Register integer data */
    OutputTime::register_data<3, FieldValue<1>::Integer>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::ELEM_DATA, &integer_field);

    /* There should be at least one output stream */
    ASSERT_EQ(OutputTime::output_streams.size(), 1);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream"));

    /* There should be two items in vector of registered element data */
    ASSERT_EQ(output_time->elem_data.size(), 2);

    /* Get first registered data */
    std::vector<OutputData*>::iterator output_data_iter = output_time->elem_data.begin();
    OutputData *output_data = *output_data_iter;

    /* There should be double data */
    ASSERT_EQ(output_data->data_type, OutputData::DOUBLE);

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->item_count, mesh.n_elements());

    /* All values has to be equal 2.0 */
    for(int i = 0; i < mesh.n_elements(); i++) {
        EXPECT_EQ(((double*)output_data->data)[i], 2.0);
    }

    /* Get next registered data */
    output_data = *(++output_data_iter);

    /* There should be double data */
    ASSERT_EQ(output_data->data_type, OutputData::INT);

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->item_count, mesh.n_elements());

    /* All values has to be equal 10 */
    for(int i = 0; i < mesh.element.size(); i++) {
        EXPECT_EQ(((int*)output_data->data)[i], 10);
    }
}
