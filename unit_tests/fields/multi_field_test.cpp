/*
 * multi_field_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>

#include <fields/multi_field.hh>
#include <fields/field_set.hh>
#include <fields/unit_si.hh>
#include <mesh/msh_gmshreader.h>
#include <input/type_base.hh>
#include <input/reader_to_storage.hh>
#include <input/type_output.hh>

#include "system/sys_profiler.hh"

#include <iostream>
using namespace std;

FLOW123D_FORCE_LINK_IN_PARENT(field_constant)

string field_constant_input = R"JSON(
{
    common={ 
        TYPE="FieldConstant", 
        value=[1, 2, 3]
    },
    transposed=[ 
        { TYPE="FieldConstant", value=1},
        { TYPE="FieldConstant", value=2},
        { TYPE="FieldConstant", value=3}
    ]
}
)JSON";

TEST(TransposeTo, field_constant) {
	MultiField<3, FieldValue<3>::Scalar> empty_mf;
	Input::Type::Record in_rec = Input::Type::Record("MultiFieldTest","")
	    .declare_key("common", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .declare_key("transposed", empty_mf.get_multifield_input_type(), Input::Type::Default::obligatory(),"" )
	    .close();

	Input::ReaderToStorage json_reader(field_constant_input, in_rec, Input::FileFormat::format_JSON);
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
	mf_assignment.set_components(component_names);
	mf_base.set_mesh( *mesh );
	mf_assignment.set_mesh( *mesh );
	mf_base.setup_components();
	mf_assignment.setup_components();

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
		EXPECT_EQ( component_names[i] + "_b", mf_assignment[i].name() );
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
		std::vector<string> component_names = { "comp_a", "comp_b", "comp_c" };
		MultiField<3, FieldValue<3>::Scalar> mf_assignment_error;
		mf_assignment_error
		    .name("d")
		    .flags(FieldFlag::input_copy);
		mf_assignment_error.set_components(component_names);
		mf_assignment_error.set_mesh( *mesh );
		mf_assignment_error.setup_components();

		EXPECT_ASSERT_DEATH( { mf_assignment_error = mf_base; },
				"Assignment between fields with different vectors of components");
	}
}
