/*
 * filed_set_test_.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: jb
 */

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include "fields/field_set.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"


enum {
	r_first,
	r_second,
	r_third
};

auto reaction_type_sel = Input::Type::Selection("ReactionType")
		.add_value(r_first, "r_first")
		.add_value(r_second, "r_second")
		.add_value(r_third, "r_third");


class SomeEquation : public testing::Test {

public:
	class EqData : public FieldSet {
	public:
		EqData() {
			*this += velocity
						.name("velocity")
						.description("Velocity vector.")
						.input_default("0.0")
						.units("LM^{-1}");
			*this += init_pressure
						.disable_where(type, {r_first, r_second })
						.name("init_pressure")
						.description("Pressure head")
						.units("L");

			*this += type
						.name("reaction_type")
						.description("")
						.input_selection(&reaction_type_sel);
		}

		// fields
	    Field<3, FieldValue<3>::Vector > velocity;
	    Field<3, FieldValue<3>::Scalar > init_pressure;
	    Field<3, FieldValue<3>::Enum > type;
	};

	SomeEquation() {
	    Profiler::initialize();
	}

	~SomeEquation() {

	}

};


TEST_F(SomeEquation, access) {
	auto data = EqData();

	EXPECT_EQ(3, data.size() );

	EXPECT_EQ(&data.velocity,  &(data["velocity"]) );
	EXPECT_EQ("velocity",  data["velocity"].name() );
	EXPECT_EQ("Velocity vector.",  data["velocity"].description() );

	FieldSet set=data.subset({"velocity", "init_pressure"});

	EXPECT_EQ(2, set.size());
	EXPECT_EQ(&data.velocity,  &(data["velocity"]) );

	EXPECT_THROW({data["noname"];}, FieldSet::ExcUnknownField);
	EXPECT_THROW({set["reaction_type"];}, FieldSet::ExcUnknownField);
}


TEST_F(SomeEquation, field_descriptor) {
	Input::Type::Record descriptor = EqData().make_field_descriptor_type("SomeEquation");

	descriptor.finish();
	EXPECT_TRUE( descriptor.has_key("time"));
	EXPECT_TRUE( descriptor.has_key("velocity"));
	EXPECT_TRUE( descriptor.has_key("init_pressure"));
	EXPECT_TRUE( descriptor.has_key("reaction_type"));
}






