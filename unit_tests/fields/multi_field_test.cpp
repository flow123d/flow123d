/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <fields/multi_field.hh>
#include <fields/field_elementwise.hh>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>

#include <mesh/msh_gmshreader.h>
#include <system/sys_profiler.hh>

#include <iostream>
using namespace std;

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

string field_constant_input = R"YAML(
common: !FieldConstant 
  value:
   - 1
   - 2
   - 3
transposed:
 - !FieldConstant
   value: 1
 - !FieldConstant
   value: 2
 - !FieldConstant
   value: 3
)YAML";

TEST(MultiField, transposition) {
	MultiField<3, FieldValue<3>::Scalar> empty_mf;
	Input::Type::Record in_rec = Input::Type::Record("MultiFieldTest","")
	    .declare_key("common", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .declare_key("transposed", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .close();

	Input::ReaderToStorage json_reader(field_constant_input, in_rec, Input::FileFormat::format_YAML);
	Input::Record input = json_reader.get_root_interface<Input::Record>();

	Input::Array common = input.val<Input::Array>("common");
	Input::Array transposed = input.val<Input::Array>("transposed");

	EXPECT_EQ(common.size(), transposed.size());

	auto it_c = common.begin<Input::AbstractRecord>();
	for (auto it_t = transposed.begin<Input::AbstractRecord>(); it_t != transposed.end(); ++it_t, ++it_c) {
		Input::Record rec_t = (*it_t);
		Input::Record rec_c = (*it_c);
		EXPECT_DOUBLE_EQ( rec_t.val<double>("value"), rec_c.val<double>("value") );
	}
}

string all_fields_input = R"YAML(
const_field: !FieldConstant
  value:
   - 1
   - 2
   - 3
elementwise_field: !FieldElementwise
  gmsh_file: fields/simplest_cube_data.msh
  field_name:
   - vector_fixed
   - vector_fixed
   - vector_fixed
formula_field: !FieldFormula
  value:
   - x
   - y-t
   - t
interpolated_p0_field: !FieldInterpolatedP0
  gmsh_file: fields/simplest_cube_data.msh
  field_name: vector_fixed
)YAML";

typedef MultiField<3, FieldValue<3>::Scalar> ScalarMultiField;
typedef ScalarMultiField::SubFieldBaseType ScalarField;

TEST(MultiField, complete_test) {
	Profiler::initialize();

	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	ScalarMultiField empty_mf;
	ScalarMultiField elementwise_mf;

	Input::Type::Record full_rec = Input::Type::Record("MultiField", "Complete multi field")
		.declare_key("const_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("elementwise_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("formula_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("interpolated_p0_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .close();

	FilePath mesh_file( "fields/simplest_cube_data.msh", FilePath::input_file);
    GmshMeshReader reader(mesh_file);
    Mesh * mesh = new Mesh;
    reader.read_physical_names(mesh);
    reader.read_mesh(mesh);
    mesh->check_and_finish();

	Input::ReaderToStorage json_reader(all_fields_input, full_rec, Input::FileFormat::format_YAML);
	Input::Record input = json_reader.get_root_interface<Input::Record>();

	Space<3>::Point point;
    point(0)=1.0; point(1)=2.0; point(2)=3.0;

    { // test of FieldConstant
		Input::Array const_fields = input.val<Input::Array>("const_field");
		EXPECT_EQ(3, const_fields.size());

		ElementAccessor<3> elm;
		double expected_val = 1.0;

		for (auto it = const_fields.begin<Input::AbstractRecord>(); it != const_fields.end(); ++it) {
			auto subfield = ScalarField::function_factory((*it), 3);
			subfield->set_time(0.0);
			auto result = subfield->value( point, elm);
			EXPECT_DOUBLE_EQ( expected_val, result );
			expected_val += 1.0;
		}
	}

	{ // test of FieldFormula
		Input::Array formula_field = input.val<Input::Array>("formula_field");
		EXPECT_EQ(3, formula_field.size());

	    ElementAccessor<3> elm;

		for (auto it = formula_field.begin<Input::AbstractRecord>(); it != formula_field.end(); ++it) {
			auto subfield = ScalarField::function_factory((*it), 3);
	        subfield->set_time(1.0);
	        auto result = subfield->value( point, elm);
			EXPECT_DOUBLE_EQ( 1.0, result );
		}
	}

	{ // test of FieldElementwise
		Input::Array elementwise_field = input.val<Input::Array>("elementwise_field");
		EXPECT_EQ(3, elementwise_field.size());

        ElementAccessor<3> el_2d=mesh->element_accessor(1);

        for (auto it = elementwise_field.begin<Input::AbstractRecord>(); it != elementwise_field.end(); ++it) {
			auto subfield = ScalarField::function_factory((*it), 3);
			subfield->set_mesh(mesh,false);
			subfield->set_time(0.0);
			auto result = subfield->value( point, el_2d);
			EXPECT_DOUBLE_EQ( 1.0, result );
		}
	}
}


