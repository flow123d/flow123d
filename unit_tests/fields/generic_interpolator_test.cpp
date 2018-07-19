/*
 * field_elementwise_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>


#include "fields/field_fe.hh"
#include "fields/generic_interpolator.hh"
#include "tools/unit_si.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "io/reader_cache.hh"

string input_str = R"YAML(
scalar: !FieldFE
  mesh_data_file: fields/interp_big_cube.msh
  field_name: scalar
  default_value: 0.0
vector_fixed: !FieldFE
  mesh_data_file: fields/interp_big_cube.msh
  field_name: vector_fixed
  default_value: 0.0
  input_discretization: native_data
)YAML";

typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
typedef FieldFE<3, FieldValue<3>::VectorFixed > VecFixField;
typedef GenericInterpolator<3, FieldValue<3>::Scalar > ScalarInterpolator;
typedef GenericInterpolator<3, FieldValue<3>::VectorFixed > VectorInterpolator;

Input::Record read_root_rec() {

    Input::Type::Record rec_type = Input::Type::Record("Test","")
        .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
        .declare_key("vector_fixed", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
        .close();

    Input::ReaderToStorage reader( input_str, rec_type, Input::FileFormat::format_YAML );
    return reader.get_root_interface<Input::Record>();
}


TEST(GenericInterpolatorTest, scalar_data) {
    std::vector<double> expected_vals = { 2.41195, 0.80715, 0.80715, 2.3756, 3.996425, 3.8303 };

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh *source_mesh = mesh_full_constructor("{mesh_file=\"fields/interp_big_cube.msh\"}");
    Mesh *target_mesh = mesh_full_constructor("{mesh_file=\"fields/interp_small_cube.msh\"}");
    Input::Record rec = read_root_rec();

    ScalarField field_in;
    ScalarField field_out;
    Space<3>::Point point;
    {
        FieldAlgoBaseInitData init_data("scalar", 0, UnitSI::dimensionless());
        field_in.init_from_input(rec.val<Input::Record>("scalar"), init_data);
        field_in.set_mesh(source_mesh,false);
        field_in.set_time(0.0);
    }
    {
        FieldAlgoBaseInitData init_data("scalar", 0, UnitSI::dimensionless());
        field_out.init_from_input(rec.val<Input::Record>("scalar"), init_data);
        field_out.set_mesh(target_mesh,false);
    }

    ScalarInterpolator interpolator;
    interpolator.interpolate(field_out, field_in);

    for(unsigned int i=0; i<target_mesh->n_elements(); i++) {
    	EXPECT_DOUBLE_EQ( expected_vals[i], field_out.value(point,target_mesh->element_accessor(i)) );
	}

}

// We need support of vector FieldFE for fix this test!
/*TEST(GenericInterpolatorTest, vector_data) {
    //std::vector<double> expected_vals = { 2.41195, 0.80715, 0.80715, 2.3756, 3.996425, 3.8303 };

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    Mesh *source_mesh = mesh_full_constructor("{mesh_file=\"fields/interp_big_cube.msh\"}");
    Mesh *target_mesh = mesh_full_constructor("{mesh_file=\"fields/interp_small_cube.msh\"}");
    Input::Record rec = read_root_rec();

    VecFixField field_in;
    VecFixField field_out;
    Space<3>::Point point;
    {
    	FieldAlgoBaseInitData init_data("vector_fixed", 0, UnitSI::dimensionless());
        field_in.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data);
        field_in.set_mesh(source_mesh,false);
        field_in.set_time(0.0);

        for(unsigned int i=0; i<source_mesh->n_elements(); i++) {
           	std::cout << " " << field_in.value(point,source_mesh->element_accessor(i));
        }
        std::cout << std::endl;
    }
    {
        FieldAlgoBaseInitData init_data("vector_fixed", 0, UnitSI::dimensionless());
        field_out.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data);
        field_out.set_mesh(target_mesh,false);
    }

    VectorInterpolator interpolator;
    interpolator.interpolate(field_out, field_in);

    for(unsigned int i=0; i<target_mesh->n_elements(); i++) {
    	EXPECT_DOUBLE_EQ( expected_vals[i], field_out.value(point,target_mesh->element_accessor(i)) );
	}

} // */
