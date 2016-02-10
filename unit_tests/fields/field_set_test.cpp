/*
 * filed_set_test_.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: jb
 */

#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "fields/field_set.hh"
#include "fields/unit_si.hh"
#include "fields/bc_field.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_constant)


enum {
	r_first,
	r_second,
	r_third
};

auto reaction_type_sel = Input::Type::Selection("ReactionType")
		.add_value(r_first, "r_first")
		.add_value(r_second, "r_second")
		.add_value(r_third, "r_third")
		.close();

// Test input for 'values' test
const string eq_data_input = R"JSON(
[
      { 
        time=0.0,
        region="BULK",
        init_pressure=1.1,
        velocity={TYPE="FieldFormula",
            value=[ "x", "y" ]
        },
        reaction_type="r_first"
      },
      { 
        time=1.0,
        region="BULK",
        velocity=[1,2]
      }
]
)JSON";

class SomeEquation : public testing::Test {

public:
	class EqData : public FieldSet {
	public:
		EqData() {
			*this += velocity
						.name("velocity")
						.description("Velocity vector.")
						.input_default("0.0")
						.flags_add(in_main_matrix)
						.units( UnitSI().kg(3).m() );
			*this += init_pressure
						.disable_where(type, {r_first, r_second })
						.name("init_pressure")
						.description("Pressure head")
						.units( UnitSI().m() );

			*this += type
						.name("reaction_type")
						.description("")
						.units( UnitSI::dimensionless() )
						.flags_add(in_main_matrix)
						.input_selection(&reaction_type_sel);
		}

		// fields
	    Field<3, FieldValue<3>::Vector > velocity;
	    Field<3, FieldValue<3>::Scalar > init_pressure;
	    Field<3, FieldValue<3>::Enum > type;
	};

	SomeEquation() {
	    Profiler::initialize();

        FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
        mesh_ = new Mesh();
        ifstream in(string(mesh_file).c_str());
        mesh_->read_gmsh_from_stream(in);
	}

	~SomeEquation() {
	    delete mesh_;
	}

	Mesh    *mesh_;
	std::vector<string> component_names_;
};


TEST_F(SomeEquation, add_operator_death) {
    auto data = EqData();
    Field<3, FieldValue<3>::Scalar > pressure;
    EXPECT_ASSERT_DEATH({
        data+=pressure
             .name("init_pressure");},
          "Another field of the same name exists");

}


TEST_F(SomeEquation, access) {
	auto data = EqData();

	// test basic add in EqData constructor
	// test size method
	// test access operator []
	EXPECT_EQ(3, data.size() );

	EXPECT_EQ(&data.velocity,  &(data["velocity"]) );
    EXPECT_EQ("velocity",  data["velocity"].name() );
    EXPECT_EQ("Velocity vector.",  data["velocity"].description() );
	EXPECT_EQ(&data.init_pressure,  &(data["init_pressure"]) );
	EXPECT_EQ(&data.type,  &(data["reaction_type"]) );

	// test subset
	FieldSet set=data.subset({"velocity", "init_pressure"});

	EXPECT_EQ(2, set.size());
	EXPECT_EQ(&data.velocity,  &(set["velocity"]) );
    EXPECT_EQ(&data.init_pressure,  &(set["init_pressure"]) );
    EXPECT_EQ(nullptr, set.field("reaction_type"));

    // test add of other field set
    Field<3, FieldValue<3>::Scalar > pressure;
    Field<3, FieldValue<3>::Scalar > copy_pressure;
    FieldSet other;

    other+=pressure.name("pressure");
    other+=copy_pressure.name("copy_pressure");
    data+=other;
    EXPECT_EQ(5, data.size());
    EXPECT_EQ(&pressure,  &(data["pressure"]) );
    EXPECT_EQ(&copy_pressure,  &(data["copy_pressure"]) );
}

TEST_F(SomeEquation, access_death) {
    auto data = EqData();

    EXPECT_THROW({data["noname"];}, FieldSet::ExcUnknownField);

}

TEST_F(SomeEquation, subset_mask) {
    auto data = EqData();

    auto main_matrix_set=data.subset(data.in_main_matrix);
    EXPECT_EQ(2, main_matrix_set.size());
    EXPECT_EQ(&data.velocity,  &(main_matrix_set["velocity"]) );
    EXPECT_EQ(&data.type,  &(main_matrix_set["reaction_type"]) );
    EXPECT_EQ(nullptr, main_matrix_set.field("init_pressure") );
}


TEST_F(SomeEquation, field_descriptor) {
	Input::Type::Record descriptor = EqData().make_field_descriptor_type("SomeEquation");

	descriptor.finish();
	EXPECT_EQ(6, descriptor.size());
	EXPECT_TRUE( descriptor.has_key("time"));
	EXPECT_TRUE( descriptor.has_key("rid"));
	EXPECT_TRUE( descriptor.has_key("region"));
	EXPECT_TRUE( descriptor.has_key("velocity"));
	EXPECT_TRUE( descriptor.has_key("init_pressure"));
	EXPECT_TRUE( descriptor.has_key("reaction_type"));
}



TEST_F(SomeEquation, output_field_selection) {
    auto data = EqData();
    BCField<3, FieldValue<3>::Scalar > bc_pressure;
    data+=bc_pressure
        .name("bc_pressure");

    Input::Type::Selection sel
        = data.make_output_field_selection("Sel", "desc").close();
    sel.finish();

    // Selection should not contain BC field bc_pressure.
    EXPECT_EQ(3, sel.size());
    EXPECT_TRUE( sel.has_name("velocity"));
    EXPECT_TRUE( sel.has_name("init_pressure"));
    EXPECT_TRUE( sel.has_name("reaction_type"));
}



TEST_F(SomeEquation, set_field) {
    auto data = EqData();
    Field<3, FieldValue<3>::Scalar > other_pressure;
    other_pressure
    .name("other_pressure")
    .description("other pressure");
    other_pressure.set_mesh(*mesh_);

    EXPECT_EQ("init_pressure", data["init_pressure"].name());
    EXPECT_EQ("Pressure head", data["init_pressure"].description());
    data["init_pressure"].flags(FieldFlag::input_copy);
    data.set_field("init_pressure", other_pressure);
    EXPECT_EQ(nullptr, data.field("other_pressure"));
    EXPECT_EQ("init_pressure", data["init_pressure"].name());
    EXPECT_EQ("other pressure", data["init_pressure"].description());

}


TEST_F(SomeEquation, collective_interface) {
    auto data = EqData();
    component_names_ = { "component_0", "component_1", "component_2", "component_3" };

    EXPECT_EQ(1,data["velocity"].n_comp());
    data.set_components(component_names_);
    EXPECT_EQ(0,data["init_pressure"].n_comp());
    EXPECT_EQ(4,data["velocity"].n_comp());
    EXPECT_EQ(0,data["reaction_type"].n_comp());

    EXPECT_EQ(nullptr,data["init_pressure"].mesh());
    data.set_mesh(*mesh_);
    EXPECT_EQ(mesh_,data["init_pressure"].mesh());
    EXPECT_EQ(mesh_,data["velocity"].mesh());
    EXPECT_EQ(mesh_,data["reaction_type"].mesh());

    // flags_add
    FieldFlag::Flags matrix(
        FieldSet::equation_input
        & FieldSet::declare_input
        & FieldSet::allow_output
        & FieldSet::in_main_matrix
        & FieldSet::in_rhs);
    FieldFlag::Flags rhs(
    FieldSet::equation_input
      & FieldSet::declare_input
      & FieldSet::allow_output
      & FieldSet::in_rhs);

    data.flags_add(FieldSet::in_rhs);
    EXPECT_EQ( rhs, data["init_pressure"].flags() );
    EXPECT_EQ( matrix, data["velocity"].flags() );
    EXPECT_EQ( matrix, data["reaction_type"].flags() );

    data.output_type(OutputTime::CORNER_DATA);
    EXPECT_EQ( OutputTime::CORNER_DATA, data["init_pressure"].output_type() );
    EXPECT_EQ( OutputTime::CORNER_DATA, data["velocity"].output_type() );
    EXPECT_EQ( OutputTime::CORNER_DATA, data["reaction_type"].output_type() );
}

TEST_F(SomeEquation, input_related) {
    auto data = EqData();
    component_names_ = { "component_0", "component_1" };

    data.set_components(component_names_);
    Input::Type::Array list_type = Input::Type::Array(data.make_field_descriptor_type("SomeEquation"));
    Input::ReaderToStorage reader( eq_data_input, list_type, Input::FileFormat::format_JSON);
    Input::Array in=reader.get_root_interface<Input::Array>();
    data.set_input_list(in);
    data.set_mesh(*mesh_);
    TimeGovernor tg(0.0, 0.5);

    data.mark_input_times(tg.equation_mark_type());
    Region front_3d = mesh_->region_db().find_id(40);
    // time = 0.0
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_FALSE(data.is_constant(front_3d));
    EXPECT_TRUE(data.changed());
    EXPECT_TRUE(tg.is_current(tg.marks().type_input()));
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_TRUE(data.changed());
    tg.next_time();

    // time = 0.5
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_FALSE(data.changed());
    EXPECT_FALSE(data.is_constant(front_3d));
    EXPECT_FALSE(tg.is_current(tg.marks().type_input()));
    tg.next_time();

    // time = 1.0
    data.set_time(tg.step(), LimitSide::right);
    EXPECT_TRUE(data.changed());
    EXPECT_TRUE(data.is_constant(front_3d));
    EXPECT_TRUE(tg.is_current(tg.marks().type_input()));

}


    /*
     * set_time
     */


    /*
     * output
     */



