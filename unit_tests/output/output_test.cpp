/*
 * vtk_output_test.cpp
 *
 *      Author: Jiri Hnidek <jiri.hnidek@tul.cz>
 */

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include "io/output.h"

#include "mesh/mesh.h"

#include "input/json_to_storage.hh"
#include "input/accessors.hh"

#include "system/sys_profiler.hh"
#include "system/file_path.hh"

#include "fields/field_base.hh"
#include "fields/field_constant.hh"


// Test #1 of input for output stream
const string output_stream1 = R"JSON(
{
  file = "./test1.pvd", 
  format = {
    TYPE = "vtk", 
    variant = "ascii"
  }, 
  name = "flow_output_stream1"
}
)JSON";


// Test #2 of input for output stream
const string output_stream2 = R"JSON(
{
  file = "./test2.pvd", 
  format = {
    TYPE = "vtk", 
    variant = "ascii"
  }, 
  name = "flow_output_stream2"
}
)JSON";


// Test input for output data
const string foo_output = R"JSON(
{
  pressure_p0 = "flow_output_stream1",
  material_id = "flow_output_stream1",
  pressure_p1 = "flow_output_stream2",
  strangeness = "flow_output_stream2"
}
)JSON";


using namespace Input::Type;


/**
 * \brief Useless class used for testing of data output
 */
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
            "Name of output stream for material ID.")
    .declare_key("pressure_p1", String(),
            "Name of output stream for P1 approximation of pressure.")
    .declare_key("strangeness", String(),
            "Name of output stream for strangeness.");


/**
 * \brief Child class used for testing OutputTime
 */
class OutputTest : public testing::Test, public OutputTime {
protected:
    virtual void SetUp() {}
    virtual void TearDown() {}
};


/**
 * \brief Test of creating of OutputTime
 */
TEST( OutputTest, test_create_output_stream ) {
    Input::JSONToStorage reader_output_stream;

    /* Read input for first output stream */
    std::stringstream in_output_stream1(output_stream1);
    reader_output_stream.read_stream(in_output_stream1, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream.get_root_interface<Input::Record>());

    /* Read input for first output stream */
    std::stringstream in_output_stream2(output_stream2);
    reader_output_stream.read_stream(in_output_stream2, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream.get_root_interface<Input::Record>());

    /* Make sure that there is one OutputTime instance */
    ASSERT_EQ(OutputTime::output_streams.size(), 2);

    /* Get first OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = (*output_iter);

    /* Pointer at OutputTime could not be NULL. Otherwise other assertation
     * would crash */
    ASSERT_TRUE(output_time != NULL);

    /* The name of instance has to be equal to configuration file:
     * "variable output_stream" */
    EXPECT_EQ(*(output_time->name), "flow_output_stream1");

    /* The type of instance should be OutputVTK */
    EXPECT_EQ(output_time->file_format, OutputTime::VTK);

    /* The filename should be "./test1.pvd" */
    EXPECT_EQ(*(output_time->base_filename()), "./test1.pvd");
}


/**
 *
 */
TEST( OutputTest, find_outputstream_by_name ) {
    /* There should be at least one output stream */
    ASSERT_EQ(OutputTime::output_streams.size(), 2);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream1"));
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

    /* There should be two output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 2);

    /* Get first OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream1"));

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

    /* Try to write data to output file */
    OutputTime::write_all_data();
}


TEST( OutputTest, test_register_corner_fields_data ) {
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
    scalar_field.set_default( Input::Type::Default("20") );
    scalar_field.set_name("pressure_p1");
    scalar_field.set_units("L");
    scalar_field.set_mesh(&mesh);
    scalar_field.set_time(0.0);

    /* Register scalar (double) data */
    OutputTime::register_data<3, FieldValue<1>::Scalar>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::CORNER_DATA, &scalar_field);

    Field<3, FieldValue<1>::Integer> integer_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    integer_field.set_default( Input::Type::Default("-1") );
    integer_field.set_name("strangeness");
    integer_field.set_units("");
    integer_field.set_mesh(&mesh);
    integer_field.set_time(0.0);

    /* Register integer data */
    OutputTime::register_data<3, FieldValue<1>::Integer>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::CORNER_DATA, &integer_field);

    /* There should two output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 2);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    /* Get second output stream */
    output_time = *(++output_iter);

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream2"));

    /* There should be two items in vector of registered element data */
    ASSERT_EQ(output_time->corner_data.size(), 2);

    /* Get first registered corner data */
    std::vector<OutputData*>::iterator output_data_iter = output_time->corner_data.begin();
    OutputData *output_data = *output_data_iter;

    /* There should be double data */
    ASSERT_EQ(output_data->data_type, OutputData::DOUBLE);

    /* All values has to be equal 20.0 */
    ElementFullIter ele = ELEMENT_FULL_ITER(&mesh, NULL);
    Node *node;
    int node_id;
    int corner_data_count, corner_id = 0;
    FOR_ELEMENTS(&mesh, ele) {
        FOR_ELEMENT_NODES(ele, node_id) {
            EXPECT_EQ(((double*)output_data->data)[corner_id], 20.0);
            corner_id++;
        }
    }

    corner_data_count = corner_id;

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->item_count, corner_data_count);

    /* Get next registered corner data */
    output_data = *(++output_data_iter);

    /* There should be double data */
    ASSERT_EQ(output_data->data_type, OutputData::INT);

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->item_count, corner_data_count);

    /* All values has to be equal 100 */
    corner_id = 0;
    FOR_ELEMENTS(&mesh, ele) {
        FOR_ELEMENT_NODES(ele, node_id) {
            EXPECT_EQ(((int*)output_data->data)[corner_id], -1);
            corner_id++;
        }
    }

    /* Try to write data to output file */
    OutputTime::write_all_data();

    OutputTime::destroy_all();
}

