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
		name = "set_1",
		region_ids= [ 2, 1 ],
		region_labels = 
		[
		   "label_3",
		   "label_2"
		] 
	},
	{
		name = "set_2",
		region_ids= 0 
	},
	{
		name = "set_3",
		union= 
		[
		   "set_1",
		   "set_2"
		] 
	}
]
)JSON";

TEST(Region, read_sets_from_input) {

	Input::JSONToStorage json_reader;
	stringstream ss(read_sets_json.c_str());
	Input::Type::Array region_set_array_input_type( RegionDB::region_set_input_type );
	json_reader.read_stream( ss,  region_set_array_input_type);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	RegionDB region_db;
	Region r0=region_db.add_region(0, "label_0", 1, false);
	Region r1=region_db.add_region(1, "label_1", 1, false);
	Region r2=region_db.add_region(2, "label_2", 2, false);
	Region r3=region_db.add_region(3, "label_3", 2, true);
	Region r4=region_db.add_region(4, "label_4", 3, false);

	region_db.read_sets_from_input(i_arr);

	EXPECT_EQ( 2,region_db.get_region_set("set_1")[0].id() );
    EXPECT_EQ( 1,region_db.get_region_set("set_1")[1].id() );
    EXPECT_EQ( 3,region_db.get_region_set("set_1")[2].id() );

    EXPECT_EQ(0, region_db.get_region_set("set_2")[0].id() );

    EXPECT_EQ(0, region_db.get_region_set("set_3")[0].id() );
    EXPECT_EQ(1, region_db.get_region_set("set_3")[1].id() );
    EXPECT_EQ(2, region_db.get_region_set("set_3")[2].id() );
    EXPECT_EQ(3, region_db.get_region_set("set_3")[3].id() );
}

