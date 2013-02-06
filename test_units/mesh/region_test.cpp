/*
 * material_dispatch_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include <gtest/gtest.h>
#include "mesh/region.hh"
#include "input/type_base.hh"
#include "input/type_output.hh"
#include "input/json_to_storage.hh"

TEST(Region, all) {
    RegionDB    region_db;

    {
    Region r=region_db.add_region(1001,"top", 2,true);
    EXPECT_EQ(2, r.idx() );
    EXPECT_TRUE( r.is_boundary() );
    EXPECT_EQ(1, r.boundary_idx() );
    EXPECT_TRUE( r.is_valid() );
    EXPECT_EQ(1001, r.id());
    EXPECT_EQ("top", r.label());
    EXPECT_EQ(2, r.dim());

    EXPECT_EQ(2, region_db.add_region(1001,"top", 2,true).idx() );
    }

    {
    Region r=region_db.add_region(1002,"inside 1", 3, false);
    EXPECT_EQ(3, r.idx() );
    EXPECT_EQ(1, r.bulk_idx() );
    EXPECT_EQ("inside 1", r.label() );
    EXPECT_EQ(1002, r.id() );
    EXPECT_EQ(3, r.dim());
    }

    {
        Region a=region_db.find_label("top");
        Region b=region_db.find_id(1001);
        EXPECT_EQ(a,b);
        Region c=region_db.find_id(1002);
        EXPECT_TRUE(a!=c);

        EXPECT_FALSE( region_db.find_id(1007).is_valid() );
    }

    {
    Region r=region_db.add_region(1003, 3);
    EXPECT_EQ(5, r.idx() );
    EXPECT_EQ(2, r.bulk_idx() );
    EXPECT_EQ("region_1003", r.label() );
    }

    {
    Region r=region_db.add_region(1004,"bottom", 2, true);
    EXPECT_EQ(4, r.idx() );
    EXPECT_EQ(2, r.boundary_idx() );
    EXPECT_EQ("bottom", r.label() );
    EXPECT_EQ(1004, r.id() );
    }

    region_db.add_region(1005,"side", 2, true);
    EXPECT_THROW( region_db.add_region(1005,"new", 3, false) , RegionDB::ExcInconsistentAdd);
    EXPECT_THROW( region_db.add_region(1001,"bottom", 2, false) , RegionDB::ExcInconsistentAdd);

    region_db.close(); // close should be called automatically at first call to any size method.

    EXPECT_EQ(8, region_db.size());
    EXPECT_EQ(4, region_db.boundary_size());
    EXPECT_EQ(3, region_db.bulk_size());

    EXPECT_DEATH( { region_db.add_region(1006,"side_", 2, true);}, "Can not add to closed region DB.");

}

const string read_sets_json = R"JSON(
[
	{
		name = "test_set_1",
		region_ids= 
		[
		   0,
		   3,
		   4
		],
		region_labels = 
		[
		   "label_1",
		   "label_2"
		] 
	},
	{
		name = "test_set_2",
		region_ids= 
		[
		   5,
		   8
		] 
	},
	{
		name = "test_set_3",
		union= 
		[
		   "set_1",
		   "set_2"
		] 
	}
]
)JSON";

TEST(Region, read_sets_from_input) {
	using namespace Input::Type;

	Record region_set_input_type =
		Record("RegionSet", "Definition of one region set.")
		.declare_key("name", String(), Default::obligatory(),
				"Unique name of the region set.")
		.declare_key("region_ids", Array( Integer(0)),
				"List of region ID numbers that has to be added to the region set.")
		.declare_key("region_labels", Array( String()),
				"List of labels of the regions that has to be added to the region set.")
		.declare_key("union", Array( String(), 2,2),
				"Defines region set as a union of given pair of sets. Overrides previous keys.")
		.declare_key("intersection", Array( String(), 2,2),
				"Defines region set as an intersection of given pair of sets. Overrides previous keys.")
		.declare_key("difference", Array( String(), 2,2),
				"Defines region set as a difference of given pair of sets. Overrides previous keys.")
		.close();

	Array region_set_array(region_set_input_type, 0, 4294967295);

	Input::JSONToStorage json_reader;
	stringstream ss(read_sets_json.c_str());
	json_reader.read_stream(ss, region_set_array);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	RegionDB region_db;
	region_db.read_sets_from_input(i_arr);

}

