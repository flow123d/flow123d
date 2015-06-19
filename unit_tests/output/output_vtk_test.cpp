/*
 * output_vtk_test.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: jb
 */

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "mesh/mesh.h"
#include "input/json_to_storage.hh"
#include "system/sys_profiler.hh"

const string test_output_time_input = R"JSON(
{
  file = "./test1.pvd", 
  format = {
    TYPE = "vtk", 
    variant = "ascii"
  } 
}
)JSON";


class TestOutputVTK : public testing::Test, public OutputVTK {
public:
    TestOutputVTK()
    : OutputVTK(
            Input::JSONToStorage(test_output_time_input, OutputTime::get_input_type())
            .get_root_interface<Input::Record>()
      )
    {
        Profiler::initialize();
        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
        this->_mesh = new Mesh();
        ifstream in(string(mesh_file).c_str());
        this->_mesh->read_gmsh_from_stream(in);
    }

    ~TestOutputVTK()
    {
        delete this->_mesh;
    }
};


TEST_F(TestOutputVTK, write_data) {
    EXPECT_EQ("./test1.pvd", this->_base_filename);
    EXPECT_EQ("test1", this->main_output_basename_);
    EXPECT_EQ(".", this->main_output_dir_);
    this->current_step=1;
    this->write_data();
}

