/*
 * output_vtk_test.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <fstream>

#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "io/output_mesh.hh"
#include "mesh/mesh.h"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "system/logger_options.hh"

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

const string test_output_time_ascii = R"YAML(
file: ./test1.pvd
format: !vtk
  variant: ascii
)YAML";

const string test_output_time_binary = R"YAML(
file: ./test1.pvd
format: !vtk
  variant: binary
)YAML";


class TestOutputVTK : public testing::Test, public OutputVTK {
public:
    TestOutputVTK()
    : OutputVTK()
    {
        Profiler::initialize();
        LoggerOptions::get_instance().set_log_file("");

        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/fields/simplest_cube_3d.msh", FilePath::input_file);
        this->_mesh = new Mesh();
        ifstream in(string(mesh_file).c_str());
        this->_mesh->read_gmsh_from_stream(in);

        component_names = { "comp_0", "comp_1", "comp_2" };
    }

    ~TestOutputVTK()
    {
        delete this->_mesh;
        LoggerOptions::get_instance().reset();
    }

    // initialize mesh with given yaml input
    void init_mesh(string input_yaml)
    {
        auto in_rec = Input::ReaderToStorage(input_yaml, OutputTime::get_input_type(), Input::FileFormat::format_YAML)
        				.get_root_interface<Input::Record>();
        this->init_from_input("dummy_equation", *(this->_mesh), in_rec);

        // create output mesh identical to computational mesh
        this->output_mesh_ = std::make_shared<OutputMesh>(*(this->_mesh));
        this->output_mesh_->create_identical_mesh();

    }

	template <class FieldType>
	void set_field_data(string field_name, string init)
    {

		// make field init it form the init string
	    FieldType field(field_name, false); // bulk field
		field.input_default(init);
		field.set_components(component_names);

		field.set_mesh( *(this->_mesh) );
		field.set_time(TimeGovernor(0.0, 1.0).step(), LimitSide::left);
		field.units(UnitSI::one());

        // create output mesh identical to computational mesh
        this->output_mesh_ = std::make_shared<OutputMesh>( *(this->_mesh) );
        this->output_mesh_->create_identical_mesh();

        this->output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>( *(this->_mesh) );
        this->output_mesh_discont_->create_mesh(this->output_mesh_);

		this->compute_field_data(ELEM_DATA, field);
	}

	// check result
	void check_result_file(std::string result_file, std::string ref_file)
	{
	    std::ifstream  vtk_file(result_file);
	    std::stringstream str_vtk_file;
	    str_vtk_file << vtk_file.rdbuf();
	    vtk_file.close();

	    std::ifstream  vtk_file_ref(ref_file);
	    std::stringstream str_vtk_file_ref;
	    str_vtk_file_ref << vtk_file_ref.rdbuf();
	    vtk_file_ref.close();

	    EXPECT_EQ(str_vtk_file_ref.str(), str_vtk_file.str());
	}

	std::vector<string> component_names;
};


TEST_F(TestOutputVTK, write_data_ascii) {
	this->init_mesh(test_output_time_ascii);
    this->current_step=1;
    this->write_data();
    EXPECT_EQ("./test1.pvd", string(this->_base_filename));
    EXPECT_EQ("test1", this->main_output_basename_);
    EXPECT_EQ(".", this->main_output_dir_);

    check_result_file("test1/test1-000001.vtu", "test_output_vtk_ascii_ref.vtu");
}


TEST_F(TestOutputVTK, write_data_binary) {
	this->init_mesh(test_output_time_binary);
    this->current_step=0;
    set_field_data< Field<3,FieldValue<0>::Scalar> > ("scalar_field", "0.5");
    set_field_data< Field<3,FieldValue<3>::VectorFixed> > ("vector_field", "[0.5, 1.0, 1.5]");
    set_field_data< Field<3,FieldValue<3>::TensorFixed> > ("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]");
    this->write_data();

    check_result_file("test1/test1-000000.vtu", "test_output_vtk_binary_ref.vtu");
}

