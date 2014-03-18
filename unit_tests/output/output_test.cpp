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
  file = "./test2.msh",
  format = {
    TYPE = "gmsh",
    variant = "ascii"
  }, 
  name = "flow_output_stream2"
}
)JSON";

// Test #1 of input for output stream
const string output_stream3 = R"JSON(
{
  file = "./test3.pvd", 
  format = {
    TYPE = "vtk", 
    variant = "ascii"
  }, 
  name = "flow_output_stream3"
}
)JSON";

// Test input for output data
const string foo_output = R"JSON(
{
  pressure_p0 = "flow_output_stream1",
  material_id = "flow_output_stream1",
  pressure_p1 = "flow_output_stream2",
  strangeness = "flow_output_stream2"
  pressure_p2 = "flow_output_stream3",
  computenode = "flow_output_stream3"
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
            "Name of output stream for strangeness.")
    .declare_key("pressure_p2", String(),
            "Name of output stream for P2 approximation of pressure.")
    .declare_key("computenode", String(),
            "Name of output stream for compute node ID.");

/**
 * \brief Child class used for testing OutputTime
 */
class OutputTest : public testing::Test {
protected:
    virtual void SetUp() {}
    virtual void TearDown() {}
    static void SetUpTestCase() { /* SetUp test case here. */ }
    static void TearDownTestCase() { OutputTime::destroy_all(); };
};


/**
 * \brief Test of creating of OutputTime
 */
TEST_F( OutputTest, test_create_output_stream ) {
    /* Read input for first output stream */
    Input::JSONToStorage reader_output_stream1(output_stream1, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream1.get_root_interface<Input::Record>());

    /* Read input for first output stream */
    Input::JSONToStorage reader_output_stream2(output_stream2, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream2.get_root_interface<Input::Record>());

    /* Read input for first output stream */
    Input::JSONToStorage reader_output_stream3(output_stream3, OutputTime::input_type);

    /* Create output stream as it is read during start of Flow */
    OutputTime::output_stream(reader_output_stream3.get_root_interface<Input::Record>());

    /* Make sure that there are 3 OutputTime instances */
    ASSERT_EQ(OutputTime::output_streams.size(), 3);

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
TEST_F( OutputTest, find_outputstream_by_name ) {
    /* There should be three output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 3);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream1"));
}


TEST_F( OutputTest, test_register_elem_fields_data ) {
    /* Read input for output */
    Input::JSONToStorage reader_output(foo_output, Foo::input_type);

    TimeGovernor tg(0.0, 1.0);

    Profiler::initialize();

    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Field<3, FieldValue<1>::Scalar> scalar_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    scalar_field.input_default( "2" );
    scalar_field.name("pressure_p0");
    scalar_field.units("L");
    scalar_field.set_mesh(mesh);
    scalar_field.set_limit_side(LimitSide::right);
    scalar_field.set_time(tg);

    /* Register scalar (double) data */
    OutputTime::register_data<3, FieldValue<1>::Scalar>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::ELEM_DATA, &scalar_field);

    Field<3, FieldValue<1>::Integer> integer_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    integer_field.input_default( "10" );
    integer_field.name("material_id");
    integer_field.units("");
    integer_field.set_mesh(mesh);
    integer_field.set_limit_side(LimitSide::right);
    integer_field.set_time(tg);

    /* Register integer data */
    OutputTime::register_data<3, FieldValue<1>::Integer>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::ELEM_DATA, &integer_field);

    /* There should be three output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 3);

    /* Get first OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream1"));

    /* There should be two items in vector of registered element data */
    ASSERT_EQ(output_time->elem_data.size(), 2);

    /* Get first registered data */
    std::vector<OutputDataBase*>::iterator output_data_iter = output_time->elem_data.begin();
    OutputDataBase *output_data = *output_data_iter;

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->items_count, mesh.n_elements());

    /* All values has to be equal 2.0 */
    for(int i = 0; i < mesh.n_elements(); i++) {
        EXPECT_EQ((*(OutputData<double>*)output_data)[i], 2.0);
    }

    /* Get next registered data */
    output_data = *(++output_data_iter);

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->items_count, mesh.n_elements());

    /* All values has to be equal 10 */
    for(int i = 0; i < mesh.element.size(); i++) {
        EXPECT_EQ((*(OutputData<int>*)output_data)[i], 10);
    }

    /* Try to write data to output file */
    OutputTime::write_all_data();
}


TEST_F( OutputTest, test_register_corner_fields_data ) {
    /* Read input for output */
    Input::JSONToStorage reader_output(foo_output, Foo::input_type);
    TimeGovernor tg(0.0, 1.0);

    Profiler::initialize();

    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Field<3, FieldValue<1>::Scalar> scalar_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    scalar_field.input_default( "20" );
    scalar_field.name("pressure_p1");
    scalar_field.units("L");
    scalar_field.set_mesh(mesh);
    scalar_field.set_limit_side(LimitSide::right);
    scalar_field.set_time(tg);

    /* Register scalar (double) data */
    OutputTime::register_data<3, FieldValue<1>::Scalar>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::CORNER_DATA, &scalar_field);

    Field<3, FieldValue<1>::Integer> integer_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    integer_field.input_default( "-1" );
    integer_field.name("strangeness");
    integer_field.units("");
    integer_field.set_mesh(mesh);
    integer_field.set_limit_side(LimitSide::right);
    integer_field.set_time(tg);

    /* Register integer data */
    OutputTime::register_data<3, FieldValue<1>::Integer>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::CORNER_DATA, &integer_field);

    /* There should three output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 3);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    /* Get second output stream */
    output_time = *(++output_iter);

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream2"));

    /* There should be two items in vector of registered element data */
    ASSERT_EQ(output_time->corner_data.size(), 2);

    /* Get first registered corner data */
    std::vector<OutputDataBase*>::iterator output_data_iter = output_time->corner_data.begin();
    OutputDataBase *output_data = *output_data_iter;

    /* All values has to be equal 20.0 */
    ElementFullIter ele = ELEMENT_FULL_ITER(&mesh, NULL);
    Node *node;
    int node_id;
    int corner_data_count, corner_id = 0;
    FOR_ELEMENTS(&mesh, ele) {
        FOR_ELEMENT_NODES(ele, node_id) {
            EXPECT_EQ((*(OutputData<double>*)output_data)[corner_id], 20.0);
            corner_id++;
        }
    }

    corner_data_count = corner_id;

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->items_count, corner_data_count);

    /* Get next registered corner data */
    output_data = *(++output_data_iter);

    /* Number of items should be equal to number of mesh elements */
    EXPECT_EQ(output_data->items_count, corner_data_count);

    /* All values has to be equal 100 */
    corner_id = 0;
    FOR_ELEMENTS(&mesh, ele) {
        FOR_ELEMENT_NODES(ele, node_id) {
            EXPECT_EQ((*(OutputData<int>*)output_data)[corner_id], -1);
            corner_id++;
        }
    }

    /* Try to write data to output file */
    OutputTime::write_all_data();

}

TEST_F( OutputTest, test_register_node_fields_data ) {
    /* Read input for output */
    Input::JSONToStorage reader_output(foo_output, Foo::input_type);
    TimeGovernor tg(0.0, 1.0);

    Profiler::initialize();

    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Field<3, FieldValue<1>::Scalar> scalar_field;

    /* Initialization of scalar field  with constant double values (1.0) */
    scalar_field.input_default( "20" );
    scalar_field.name("pressure_p2");
    scalar_field.units("L");
    scalar_field.set_mesh(mesh);
    scalar_field.set_limit_side(LimitSide::right);
    scalar_field.set_time(tg);

    /* Register scalar (double) data */
    OutputTime::register_data<3, FieldValue<1>::Scalar>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::NODE_DATA, &scalar_field);

    Field<3, FieldValue<1>::Integer> integer_field;

    /* Initialization of scalar field  with constant int values (0) */
    integer_field.input_default( "2" );
    integer_field.name("computenode");
    integer_field.units("");
    integer_field.set_mesh(mesh);
    integer_field.set_limit_side(LimitSide::right);
    integer_field.set_time(tg);

    /* Register integer data */
    OutputTime::register_data<3, FieldValue<1>::Integer>(reader_output.get_root_interface<Input::Record>(),
            OutputTime::NODE_DATA, &integer_field);

    /* There should three output streams */
    ASSERT_EQ(OutputTime::output_streams.size(), 3);

    /* Get this OutputTime instance */
    std::vector<OutputTime*>::iterator output_iter = OutputTime::output_streams.begin();
    OutputTime *output_time = *output_iter;

    /* Get third output stream */
    ++output_iter;
    output_time = *(++output_iter);

    ASSERT_EQ(output_time, OutputTime::output_stream_by_name("flow_output_stream3"));

    /* There should be two items in vector of registered element data */
    ASSERT_EQ(output_time->node_data.size(), 2);

    /* Get first registered corner data */
    std::vector<OutputDataBase*>::iterator output_data_iter = output_time->node_data.begin();
    OutputDataBase *output_data = *output_data_iter;

    /* Number of items should be equal to number of mesh nodes */
    EXPECT_EQ(output_data->items_count, mesh.n_nodes());

    /* All values has to be equal 20.0 */
    for(int node_id=0; node_id < mesh.n_nodes(); node_id++) {
        EXPECT_EQ((*(OutputData<double>*)output_data)[node_id], 20.0);
    }

    /* Get next registered corner data */
    output_data = *(++output_data_iter);

    /* Number of items should be equal to number of mesh nodes */
    EXPECT_EQ(output_data->items_count, mesh.n_nodes());

    /* All values has to be equal 2 */
    for(int node_id=0; node_id < mesh.n_nodes(); node_id++) {
        EXPECT_EQ((*(OutputData<int>*)output_data)[node_id], 2);
    }

    /* Try to write data to output file */
    OutputTime::write_all_data();

}

