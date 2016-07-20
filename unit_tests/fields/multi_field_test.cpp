/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include <fields/multi_field.hh>
#include <fields/field_elementwise.hh>
#include <fields/field_set.hh>
#include <fields/unit_si.hh>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>
#include <input/type_output.hh>
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


string eq_data_input = R"JSON(
[
  { id=37,
    a=1,
    b=0
  }
] 
}
)JSON";

TEST(Operators, assignment) {
    Profiler::initialize();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
    GmshMeshReader msh_reader(mesh_file);
    Mesh * mesh = new Mesh;
    msh_reader.read_physical_names(mesh);
	msh_reader.read_mesh(mesh);
	mesh->check_and_finish();

	std::vector<string> component_names = { "comp_0", "comp_1", "comp_2" };

	MultiField<3, FieldValue<3>::Scalar> mf_base;
	mf_base.name("a");
	mf_base.set_components(component_names);
	mf_base.units( UnitSI::dimensionless() );

	FieldSet set_of_field;
	set_of_field += mf_base;
    Input::Type::Array list_type = Input::Type::Array(set_of_field.make_field_descriptor_type("MultiFieldTest"));
    Input::ReaderToStorage reader( eq_data_input, list_type, Input::FileFormat::format_JSON);
    Input::Array in_list=reader.get_root_interface<Input::Array>();
    set_of_field.set_input_list(in_list);

    MultiField<3, FieldValue<3>::Scalar> mf_assignment;
	EXPECT_EQ("", mf_assignment.name());
	EXPECT_FALSE(mf_assignment.is_bc());

	// copies
	mf_assignment
	    .name("b")
	    .flags(FieldFlag::input_copy);
	std::vector<string> component_names_2 = { "comp_a", "comp_b", "comp_c" };
	mf_assignment.set_components(component_names_2);
	mf_base.set_mesh( *mesh );
	mf_assignment.set_mesh( *mesh );
	mf_base.setup_components();

	MultiField<3, FieldValue<3>::Scalar> mf_copy(mf_base);	// copy constructor
	mf_assignment = mf_base; // assignment

	EXPECT_STREQ("a", mf_base.name().c_str());
	EXPECT_STREQ("b", mf_assignment.name().c_str());
	EXPECT_STREQ("a", mf_copy.name().c_str());
	EXPECT_EQ(3, mf_base.size());
	EXPECT_EQ(mf_assignment.size(), mf_base.size());
	EXPECT_EQ(mf_copy.size(), mf_base.size());
	for (unsigned int i=0; i<mf_base.size(); ++i) {
		EXPECT_EQ( component_names[i] + "_a", mf_base[i].name() );
		EXPECT_EQ( component_names_2[i] + "_b", mf_assignment[i].name() );
		EXPECT_EQ( mf_base[i].name(), mf_copy[i].name() );
	}

	{
		// throw assert, empty vector of component names
		MultiField<3, FieldValue<3>::Scalar> mf_assignment_error;
		mf_assignment_error
		    .name("c")
		    .flags(FieldFlag::input_copy);
		mf_assignment_error.set_mesh( *mesh );

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Vector of component names can't be empty!");
	}

	{
		// throw assert, source field has different component names
		std::vector<string> component_names = { "comp_a", "comp_b" };
		MultiField<3, FieldValue<3>::Scalar> mf_assignment_error;
		mf_assignment_error
		    .name("d")
		    .flags(FieldFlag::input_copy);
		mf_assignment_error.set_components(component_names);
		mf_assignment_error.set_mesh( *mesh );
		mf_assignment_error.setup_components();

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Both multi fields must have same size of vectors");
	}
}
