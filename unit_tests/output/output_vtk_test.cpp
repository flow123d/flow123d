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

class TestOutputVTK : public OutputVTK, public std::enable_shared_from_this<OutputVTK> {
public:
    TestOutputVTK()
    : OutputVTK()
    {
        LoggerOptions::get_instance().set_log_file("");

        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/fields/simplest_cube_3d.msh", FilePath::input_file);
        this->_mesh = mesh_full_constructor("{ mesh_file=\"" + (string)mesh_file + "\", optimize_mesh=false }");

        component_names = { "comp_0", "comp_1", "comp_2" };

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
		Field<spacedim, Value> field(field_name); // bulk field
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

	template <class FieldVal>
	void set_native_field_data(string field_name, unsigned int size, double step)
    {

		// make field init it form the init string
		Field<3, FieldVal> field(field_name);
		field.set_components(component_names);
		field.set_mesh( *(this->_mesh) );
		field.units(UnitSI::one());

		std::shared_ptr<DOFHandlerMultiDim> dh = make_shared<DOFHandlerMultiDim>( *(this->_mesh) );
		MixedPtr<FE_P_disc> fe(0);
        std::shared_ptr<::DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(this->_mesh, fe);
		dh->distribute_dofs(ds);

		VectorMPI v(size);
        for (unsigned int i=0; i<size; ++i) v.set(i, step*i);

        auto native_data_ptr = make_shared< FieldFE<3, FieldVal> >();
        native_data_ptr->set_fe_data(dh, v);

        field.set(native_data_ptr, 0.0);
        field.output_type(OutputTime::NATIVE_DATA);
        field.set_time(TimeGovernor(0.0, 1.0).step(), LimitSide::left);

        field.compute_field_data(OutputTime::NATIVE_DATA, shared_from_this());
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

	std::string base_filename() {
		return string(this->_base_filename);
	}

	std::string main_output_basename() {
		return this->main_output_basename_;
	}

	std::string main_output_dir() {
		return this->main_output_dir_;
	}

	std::vector<string> component_names;
	Mesh *_mesh;
	std::shared_ptr<OutputMeshBase> output_mesh_;
};


TEST_F(TestVTK, write_data_ascii) {
	std::shared_ptr<TestOutputVTK> output_vtk = std::make_shared<TestOutputVTK>();

	output_vtk->init_mesh(test_output_time_ascii);
	output_vtk->set_current_step(0);
	output_vtk->set_field_data<3, FieldValue<0>::Scalar> ("scalar_field", "0.5", "0.5");
	output_vtk->set_field_data<3, FieldValue<3>::VectorFixed> ("vector_field", "[0.5, 1.0, 1.5]", "0.5 1.0 1.5");
	output_vtk->set_field_data<3, FieldValue<3>::TensorFixed> ("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]", "1 2 3; 4 5 6; 7 8 9");
	output_vtk->set_native_field_data< FieldValue<0>::Scalar >("flow_data", 6, 0.2);
	// output_vtk->write_data();
	// output_vtk->check_result_file("test1/test1-000000.vtu", "test_output_vtk_ascii_ref.vtu");

	output_vtk->clear_data();
	output_vtk->set_current_step(1);
	output_vtk->set_field_data<3, FieldValue<0>::Scalar> ("scalar_field", "0.5", "0.5");
	output_vtk->set_field_data<3, FieldValue<3>::VectorFixed> ("vector_field", "[0.5, 1.0, 1.5]", "0.5 1.0 1.5");
	output_vtk->set_field_data<3, FieldValue<3>::TensorFixed> ("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]", "1 2 3; 4 5 6; 7 8 9");
	output_vtk->set_native_field_data< FieldValue<0>::Scalar >("flow_data", 6, 0.2);
	output_vtk->write_data();

	EXPECT_EQ("./test1.pvd", output_vtk->base_filename());
    EXPECT_EQ("test1", output_vtk->main_output_basename());
    EXPECT_EQ(".", output_vtk->main_output_dir());
	output_vtk->check_result_file("test1/test1-000001.vtu", "test_output_vtk_ascii_ref.vtu");
}


TEST_F(TestVTK, write_data_binary) {
	std::shared_ptr<TestOutputVTK> output_vtk = std::make_shared<TestOutputVTK>();

	output_vtk->init_mesh(test_output_time_binary);
    output_vtk->set_current_step(0);
    output_vtk->set_field_data<3, FieldValue<0>::Scalar>("scalar_field", "0.5", "0.5");
    output_vtk->set_field_data<3, FieldValue<3>::VectorFixed>("vector_field", "[0.5, 1.0, 1.5]", "0.5 1.0 1.5");
    output_vtk->set_field_data<3, FieldValue<3>::TensorFixed>("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]", "1 2 3; 4 5 6; 7 8 9");
    output_vtk->write_data();

    output_vtk->check_result_file("test1/test1-000000.vtu", "test_output_vtk_binary_ref.vtu");
}

#ifdef FLOW123D_HAVE_ZLIB

const string test_output_time_compressed = R"YAML(
file: ./test1.pvd
format: !vtk
  variant: binary_zlib
)YAML";

TEST_F(TestVTK, write_data_compressed) {
	std::shared_ptr<TestOutputVTK> output_vtk = std::make_shared<TestOutputVTK>();

	output_vtk->init_mesh(test_output_time_compressed);
    output_vtk->set_current_step(0);
    output_vtk->set_field_data<3, FieldValue<0>::Scalar>("scalar_field", "0.5", "0.5");
    output_vtk->set_field_data<3, FieldValue<3>::VectorFixed>("vector_field", "[0.5, 1.0, 1.5]", "0.5 1.0 1.5");
    output_vtk->set_field_data<3, FieldValue<3>::TensorFixed>("tensor_field", "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]", "1 2 3; 4 5 6; 7 8 9");
    output_vtk->write_data();

    output_vtk->check_result_file("test1/test1-000000.vtu", "test_output_vtk_zlib_ref.vtu");
}

#endif // FLOW123D_HAVE_ZLIB

