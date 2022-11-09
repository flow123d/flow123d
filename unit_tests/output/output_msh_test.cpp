/*
 * output_vtk_test.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <fstream>

#include "config.h"

#include "io/output_time.hh"
#include "io/output_msh.hh"
#include "io/output_mesh.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "input/reader_to_storage.hh"
#include "system/logger_options.hh"
#include "system/sys_profiler.hh"
#include "fields/field.hh"

const string test_output_time = R"YAML(
file: ./test_output.msh
format: !gmsh
  variant: ascii
)YAML";


class TestOutputMSH : public OutputMSH, public std::enable_shared_from_this<OutputMSH> {
public:
	TestOutputMSH()
    : OutputMSH()
    {
        Profiler::instance();
        LoggerOptions::get_instance().set_log_file("");

        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/fields/simplest_cube_3d.msh", FilePath::input_file);
        this->_mesh = mesh_full_constructor("{ mesh_file=\"" + (string)mesh_file + "\", optimize_mesh=false }");

        component_names = { "comp_0", "comp_1", "comp_2" };

        this->write_time = 0.0; // hack: unset condition in OutputTime::write_time_frame and output is not performed
    }

    ~TestOutputMSH()
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
        output_mesh_ = std::make_shared<OutputMesh>(*(this->_mesh));
        output_mesh_->create_sub_mesh();
        output_mesh_->make_serial_master_mesh();
        this->set_output_data_caches(output_mesh_);

    }

	template <int spacedim, class Value>
	void set_field_data(string field_name, string init, string rval)
    {
		typedef typename Value::element_type ElemType;

		// make field, init it form the init string
		Field<spacedim, Value> field(field_name, false); // bulk field
		field.input_default(init);
		field.set_components(component_names);

		field.set_mesh( *(this->_mesh) );
		field.units(UnitSI::one());
		field.set_time(TimeGovernor(0.0, 1.0).step(), LimitSide::left);
		field.set_output_data_cache(OutputTime::ELEM_DATA, shared_from_this());
	    auto output_cache_base = this->prepare_compute_data<ElemType>(field_name, OutputTime::ELEM_DATA,
	            (unsigned int)Value::NRows_, (unsigned int)Value::NCols_);
	    std::shared_ptr<ElementDataCache<ElemType>> output_data_cache = std::dynamic_pointer_cast<ElementDataCache<ElemType>>(output_cache_base);
	    arma::mat ret_value(rval);
	    for (uint i=0; i<output_data_cache->n_values(); ++i)
	        output_data_cache->store_value(i, ret_value.memptr() );

        //this->output_mesh_discont_ = std::make_shared<OutputMeshDiscontinuous>( *(this->_mesh) );
        //this->output_mesh_discont_->create_sub_mesh();
        //this->output_mesh_discont_->make_serial_master_mesh();

	    this->update_time(field.time());
	}

	// check result
    void check_result_file(std::string result_file, std::string ref_file)
    {
        std::ifstream  msh_file(result_file);
        std::stringstream str_msh_file;
        str_msh_file << msh_file.rdbuf();
        msh_file.close();

        std::ifstream  msh_file_ref(ref_file);
        std::stringstream str_msh_file_ref;
        str_msh_file_ref << msh_file_ref.rdbuf();
        msh_file_ref.close();

        EXPECT_EQ(str_msh_file_ref.str(), str_msh_file.str());
    }

	void set_current_step(int step) {
		this->current_step = step;
	}

	std::string base_filename() {
		return string(this->_base_filename);
	}

	std::vector<string> component_names;
	Mesh *_mesh;
	std::shared_ptr<OutputMeshBase> output_mesh_;
};


TEST(TestOutputMSH, write_data) {
	std::shared_ptr<TestOutputMSH> output_msh = std::make_shared<TestOutputMSH>();
	output_msh->init_mesh(test_output_time);

	output_msh->set_current_step(0);
	output_msh->set_field_data<3, FieldValue<0>::Scalar> ("scalar_field", "0.5", "0.5");
	output_msh->set_field_data<3, FieldValue<3>::Vector> ("vector_field", "[0.5, 1.0, 1.5]", "0.5 1.0 1.5");
	output_msh->set_field_data<3, FieldValue<3>::Tensor> ("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]", "1 2 3; 4 5 6; 7 8 9");
	output_msh->write_data();

	output_msh->clear_data();
	output_msh->set_current_step(1);
	output_msh->set_field_data<3, FieldValue<0>::Scalar> ("scalar_field", "0.75", "0.75");
	output_msh->set_field_data<3, FieldValue<3>::Vector> ("vector_field", "[0.75, 1.5, 2.25]", "0.75 1.5 2.25");
	output_msh->set_field_data<3, FieldValue<3>::Tensor> ("tensor_field", "[[1, 4, 7], [2, 5, 8], [3, 6, 9]]", "1 4 7; 2 5 8; 3 6 9");
	output_msh->write_data();

    EXPECT_EQ("./test_output.msh", output_msh->base_filename());
    output_msh->check_result_file("./test_output.msh", "./test_output_gmsh_ref.msh");
}
