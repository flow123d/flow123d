/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include <fields/multi_field.hh>
#include <fields/field_fe.hh>
#include <fields/field_set.hh>
#include <tools/unit_si.hh>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>
#include <input/type_output.hh>
#include <io/msh_gmshreader.h>
#include <system/sys_profiler.hh>

#include <iostream>
using namespace std;

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)
FLOW123D_FORCE_LINK_IN_PARENT(field_formula)
FLOW123D_FORCE_LINK_IN_PARENT(field_fe)

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
const_field_full: !FieldConstant
  value:
   - 1
   - 2
   - 3
const_field_base: !FieldConstant
  value: 1
const_field_autoconv: 1

formula_field_full: !FieldFormula
  value:
   - t
   - x
   - y-t
formula_field_base: !FieldFormula
  value: x

fe_field: !FieldFE
  mesh_data_file: fields/simplest_cube_data.msh
  field_name: vector_fixed
interpolated_p0_field: !FieldFE
  mesh_data_file: fields/simplest_cube_3d.msh
  field_name: scalar
  interpolation: P0_intersection
)YAML";

class MultiFieldTest : public testing::Test {
public:
	typedef MultiField<3, FieldValue<3>::Scalar> ScalarMultiField;
	typedef ScalarMultiField::SubFieldBaseType ScalarField;

protected:
    virtual void SetUp() {
    	Profiler::instance();
    	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    	point(0)=1.0; point(1)=2.0; point(2)=3.0;

        auto mesh_reader = reader_constructor("{mesh_file=\"fields/simplest_cube_data.msh\"}");
        mesh = mesh_constructor("{mesh_file=\"fields/simplest_cube_data.msh\"}");
        mesh_reader->read_raw_mesh(mesh);
        mesh->setup_topology();

    }

    virtual void TearDown() {
    	delete mesh;
    }

    void check_field_vals(Input::Array &arr_field, ElementAccessor<3> elm, double expected = 1.0, double step = 0.0) {
    	for (auto it = arr_field.begin<Input::AbstractRecord>(); it != arr_field.end(); ++it) {
    	    FieldAlgoBaseInitData init_data("test_mf", 3, UnitSI::dimensionless());
    		auto subfield = ScalarField::function_factory((*it), init_data);
    		subfield->set_mesh(mesh, false);
    		subfield->set_time(0.0);
    		auto result = subfield->value( point, elm );
    		EXPECT_DOUBLE_EQ( expected, result );
    		expected += step;
    	}
    }

    static Input::Type::Record & get_input_type();
    static ScalarMultiField empty_mf;
    Mesh * mesh;
    Space<3>::Point point;
};

MultiFieldTest::ScalarMultiField MultiFieldTest::empty_mf = MultiFieldTest::ScalarMultiField();

Input::Type::Record & MultiFieldTest::get_input_type() {
	return Input::Type::Record("MultiField", "Complete multi field")
		.declare_key("const_field_full", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("const_field_base", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("const_field_autoconv", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("formula_field_full", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("formula_field_base", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("fe_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.declare_key("interpolated_p0_field", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
		.close();
}


TEST_F(MultiFieldTest, complete_test) {
	Input::ReaderToStorage json_reader(all_fields_input, MultiFieldTest::get_input_type(), Input::FileFormat::format_YAML);
	Input::Record input = json_reader.get_root_interface<Input::Record>();

    { // test of FieldConstant - full input
		Input::Array const_fields = input.val<Input::Array>("const_field_full");
		EXPECT_EQ(3, const_fields.size());

		ElementAccessor<3> elm;
		check_field_vals(const_fields, elm, 1.0, 1.0);
	}

    { // test of FieldConstant - set key 'value' with one value
		Input::Array const_fields = input.val<Input::Array>("const_field_base");
		EXPECT_EQ(1, const_fields.size());

		ElementAccessor<3> elm;
		check_field_vals(const_fields, elm);
	}

    { // test of FieldConstant - autoconversion of FieldAlgorithmBase Abstract and 'value' key
		Input::Array const_fields = input.val<Input::Array>("const_field_autoconv");
		EXPECT_EQ(1, const_fields.size());

		ElementAccessor<3> elm;
		check_field_vals(const_fields, elm);
	}

	{ // test of FieldFormula - full input
		Input::Array formula_field = input.val<Input::Array>("formula_field_full");
		EXPECT_EQ(3, formula_field.size());

	    ElementAccessor<3> elm;
	    check_field_vals(formula_field, elm, 0.0, 1.0);
	}

	{ // test of FieldFormula - set key 'value' with one value
		Input::Array formula_field = input.val<Input::Array>("formula_field_base");
		EXPECT_EQ(1, formula_field.size());

	    ElementAccessor<3> elm;
	    check_field_vals(formula_field, elm);
	}

	{ // test of FieldFE
		Input::Array fe_field = input.val<Input::Array>("fe_field");
		EXPECT_EQ(1, fe_field.size());

        check_field_vals(fe_field, mesh->element_accessor(1));
	}

	{ // test of FieldInterpolatedP0
		Input::Array interpolated_p0_field = input.val<Input::Array>("interpolated_p0_field");
		EXPECT_EQ(1, interpolated_p0_field.size());

        check_field_vals(interpolated_p0_field, mesh->element_accessor(1), 0.650, 0.0);
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
    Profiler::instance();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    auto mesh_reader = reader_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    Mesh * mesh = mesh_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");
    mesh_reader->read_physical_names(mesh);
    mesh_reader->read_raw_mesh(mesh);
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
    TimeGovernor tg(0.0, 1.0);
    set_of_field.set_input_list(in_list, tg);

    MultiField<3, FieldValue<3>::Scalar> mf_assignment;
	EXPECT_EQ("", mf_assignment.name());
	EXPECT_FALSE(mf_assignment.is_bc());

	// copies
	mf_assignment
	    .name("b")
	    .flags(FieldFlag::input_copy)
	    .units( UnitSI::dimensionless() );
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
		    .flags(FieldFlag::input_copy)
		    .units( UnitSI::dimensionless() );
		mf_assignment_error.set_components(component_names);
		mf_assignment_error.set_mesh( *mesh );
		mf_assignment_error.setup_components();

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Both multi fields must have same size of vectors");
	}

	delete mesh;
}
