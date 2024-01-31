/*
 * optimizer_output_test.cpp
 *
 * Perform output of optimized element indices.
 *
 *  Created on: Dec 2, 2020
 *      Author: df
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <fstream>

#include "config.h"

#include "io/output_time.hh"
#include "io/output_vtk.hh"
#include "io/output_mesh.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "input/reader_to_storage.hh"
#include "system/logger_options.hh"
#include "system/sys_profiler.hh"
#include "fields/field.hh"

#include "fem/mapping_p1.hh"
#include "fem/dofhandler.hh"
#include "fem/fe_p.hh"
#include "fields/field_fe.hh"
#include "la/vector_mpi.hh"
#include "tools/mixed.hh"

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

const string test_output_time = R"YAML(
file: ./test1.pvd
format: !vtk
  variant: ascii
)YAML";


class TestVTK : public testing::Test {
protected:
	TestVTK()
    {
        Profiler::instance();
    }

    ~TestVTK()
    {
        Profiler::uninitialize();
    }
};

class TestOutputVTK : public OutputVTK {
public:
    TestOutputVTK()
    : OutputVTK()
    {
        LoggerOptions::get_instance().set_no_log();

        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/../tests/00_mesh/square_1x1_frac_fork.msh", FilePath::input_file);
        this->_mesh = mesh_full_constructor("{ mesh_file=\"" + (string)mesh_file + "\" }");

        this->write_time = 0.0; // hack: unset condition in OutputTime::write_time_frame and output is not performed
    }

    ~TestOutputVTK()
    {
        delete this->_mesh;
        LoggerOptions::get_instance().reset();
    }

    // initialize mesh with given yaml input
    void init_mesh(string input_yaml)
    {
    	auto in_rec = Input::ReaderToStorage(input_yaml, const_cast<Input::Type::Record &>(OutputTime::get_input_type()), Input::FileFormat::format_YAML)
        				.get_root_interface<Input::Record>();
        this->init_from_input("dummy_equation", in_rec, std::make_shared<TimeUnitConversion>());

        // create output mesh identical to computational mesh
        auto output_mesh = std::make_shared<OutputMesh>(*(this->_mesh));
        output_mesh->create_sub_mesh();
        output_mesh->make_serial_master_mesh();
        this->set_output_data_caches(output_mesh);

    }

    void set_mesh_optimized_field() {
        OutputTime::OutputDataPtr output_data_base = this->prepare_compute_data<unsigned int>("hilbert", ELEM_DATA,
        		(unsigned int)1, (unsigned int)1);

        auto output_data = std::dynamic_pointer_cast<ElementDataCache<unsigned int>>(output_data_base);
        auto &data_vec = *( output_data->get_data().get() );
        for (uint i=0; i<data_vec.size(); ++i) {
            data_vec[i] = i;
        }
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

	void set_current_step(int step) {
		this->current_step = step;
	}

	Mesh *_mesh;
	std::shared_ptr<OutputMeshBase> output_mesh_;
};


TEST_F(TestVTK, hilbert) {
    std::shared_ptr<TestOutputVTK> output_vtk = std::make_shared<TestOutputVTK>();

    output_vtk->init_mesh(test_output_time);
    output_vtk->set_current_step(0);
    output_vtk->set_mesh_optimized_field();
    output_vtk->write_data();
    // Check output './test1/test1-000000.vtu' visually in Paraview.
}

