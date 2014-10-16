/*
 * material_dispatch_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include <flow_gtest.hh>
#include "mesh/region.hh"
#include "input/type_base.hh"
#include "input/type_output.hh"
#include "input/json_to_storage.hh"
#include <map>

#include <boost/lexical_cast.hpp>

TEST(Region, all) {
//    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    RegionDB    region_db;
    region_db.add_region(0,".nothing_bc", 0);
    region_db.add_region(1,"nothing_bulk", 0);

    {
    Region r=region_db.add_region(1001,".top", 2);
    EXPECT_EQ(2, r.idx() );
    EXPECT_TRUE( r.is_boundary() );
    EXPECT_EQ(1, r.boundary_idx() );
    EXPECT_TRUE( r.is_valid() );
    EXPECT_EQ(1001, r.id());
    EXPECT_EQ(".top", r.label());
    EXPECT_EQ(2, r.dim());

    EXPECT_EQ(2, region_db.add_region(1001,".top", 2).idx() );

    // try to convert Region to RegionIdx
    RegionIdx r_idx = r;
    EXPECT_EQ(2, r_idx.idx());
    }

    {
    Region r=region_db.add_region(1002,"inside 1", 3);
    EXPECT_EQ(3, r.idx() );
    EXPECT_EQ(1, r.bulk_idx() );
    EXPECT_EQ("inside 1", r.label() );
    EXPECT_EQ(1002, r.id() );
    EXPECT_EQ(3, r.dim());
    }

    {
        Region a=region_db.find_label(".top");
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
    Region r=region_db.add_region(1004,".bottom", 2);
    EXPECT_EQ(4, r.idx() );
    EXPECT_EQ(2, r.boundary_idx() );
    EXPECT_EQ(".bottom", r.label() );
    EXPECT_EQ(1004, r.id() );
    }

    region_db.add_region(1005,".side", 2);
    //EXPECT_THROW( region_db.add_region(1005,"new", 3, false) , RegionDB::ExcInconsistentAdd);
    EXPECT_THROW( region_db.add_region(1001,".bottom", 2) , RegionDB::ExcNonuniqueID);
    EXPECT_THROW( region_db.add_region(1101,".top", 2) , RegionDB::ExcNonuniqueLabel);
    EXPECT_THROW( region_db.add_region(1001,".top", 3) , RegionDB::ExcNonuniqueLabel);
    //EXPECT_THROW( region_db.add_region(1001,"top", 2) , RegionDB::ExcNonuniqueLabel);

    region_db.close(); // close should be called automatically at first call to any size method.

    EXPECT_EQ(8, region_db.size());
    EXPECT_EQ(4, region_db.boundary_size());
    EXPECT_EQ(3, region_db.bulk_size());

    RegionSet bulk = region_db.get_region_set("BULK");
    EXPECT_EQ(3,bulk.size());
    EXPECT_EQ(1, bulk[0].id());
    EXPECT_EQ(1002, bulk[1].id());
    EXPECT_EQ(1003, bulk[2].id());

    EXPECT_THROW_WHAT( { region_db.add_region(1006,".side_", 2);}, RegionDB::ExcAddingIntoClosed, "Can not add label='.side_'");

}

TEST(Region, add_nonunique_id_region) {
	RegionDB region_db;

	{
		Region r=region_db.add_region(1, "user_defined_name");
	    EXPECT_EQ(1, r.idx() );
	    EXPECT_FALSE( r.is_boundary() );
	    EXPECT_EQ(0, r.bulk_idx() );
	    EXPECT_TRUE( r.is_valid() );
	    EXPECT_EQ(1, r.id());
	    EXPECT_EQ("user_defined_name", r.label());
	}

	{
		Region r=region_db.add_region(1, "region_name", 3);
	    EXPECT_EQ(1, r.idx() );
	    EXPECT_FALSE( r.is_boundary() );
	    EXPECT_EQ(0, r.bulk_idx() );
	    EXPECT_TRUE( r.is_valid() );
	    EXPECT_EQ(1, r.id());
	    EXPECT_EQ("user_defined_name", r.label());
	    EXPECT_EQ(3, r.dim());
	}

	{
		Region r=region_db.add_region(2, "other_dim", 2);
	    EXPECT_EQ(3, r.idx() );
	    EXPECT_FALSE( r.is_boundary() );
	    EXPECT_EQ(1, r.bulk_idx() );
	    EXPECT_TRUE( r.is_valid() );
	    EXPECT_EQ(2, r.id());
	    EXPECT_EQ("other_dim", r.label());
	    EXPECT_EQ(2, r.dim());
	}

}

const string read_sets_json = R"JSON(
[
	{
		name = "set_1",
		region_ids= [ 2, 1 ],
		region_labels = 
		[
		   ".label_3",
		   "label_2"
		] 
	},
	{
		name = "set_2",
		region_ids= 0 
	},
	{
		name = "set_3",
		region_ids= [ 4, 2, 3 ] 
	},
	{
		name = "set_4",
		union= 
		[
		   "set_1",
		   "set_2"
		] 
	},
	{
		name = "set_5",
		difference= 
		[
		   "set_1",
		   "set_3"
		] 
	},
	{
		name = "set_6",
		intersection= 
		[
		   "set_1",
		   "set_3"
		] 
	}
]
)JSON";

TEST(Region, read_sets_from_input) {
	Input::Type::Array region_set_array_input_type( RegionDB::region_set_input_type );
	Input::JSONToStorage json_reader( read_sets_json,  region_set_array_input_type);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	RegionDB region_db;
	region_db.add_region(0, "label_0", 1);
	region_db.add_region(1, "label_1", 1);
	region_db.add_region(2, "label_2", 2);
	region_db.add_region(3, ".label_3", 2);
	region_db.add_region(4, "label_4", 3);

	region_db.read_sets_from_input(i_arr);

	EXPECT_EQ(3, region_db.get_region_set("set_1").size() );
	EXPECT_EQ(2, region_db.get_region_set("set_1")[0].id() );
    EXPECT_EQ(1, region_db.get_region_set("set_1")[1].id() );
    EXPECT_EQ(3, region_db.get_region_set("set_1")[2].id() );

    EXPECT_EQ(1, region_db.get_region_set("set_2").size() );
    EXPECT_EQ(0, region_db.get_region_set("set_2")[0].id() );

	EXPECT_EQ(3, region_db.get_region_set("set_3").size() );
	EXPECT_EQ(4, region_db.get_region_set("set_3")[0].id() );
    EXPECT_EQ(2, region_db.get_region_set("set_3")[1].id() );
    EXPECT_EQ(3, region_db.get_region_set("set_3")[2].id() );

    EXPECT_EQ(4, region_db.get_region_set("set_4").size() );
    EXPECT_EQ(3, region_db.get_region_set("set_4")[0].id() );
    EXPECT_EQ(0, region_db.get_region_set("set_4")[1].id() );
    EXPECT_EQ(1, region_db.get_region_set("set_4")[2].id() );
    EXPECT_EQ(2, region_db.get_region_set("set_4")[3].id() );

    EXPECT_EQ(1, region_db.get_region_set("set_5").size() );
    EXPECT_EQ(1, region_db.get_region_set("set_5")[0].id() );

    EXPECT_EQ(2, region_db.get_region_set("set_6").size() );
    EXPECT_EQ(3, region_db.get_region_set("set_6")[0].id() );
    EXPECT_EQ(2, region_db.get_region_set("set_6")[1].id() );
}

const string read_element_map_json = R"JSON(
[
	{
		name = "label_0",
		id = 0,
		element_list = [0, 5, 9]
	},
	{
		name = "label_1",
		id = 1,
		element_list = [8, 5, 3, 1]
	}
]
)JSON";

TEST(Region, read_element_map_from_input) {

	Input::Type::Array element_map_array_input_type( RegionDB::region_input_type );
	Input::JSONToStorage json_reader( read_element_map_json,  element_map_array_input_type);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	RegionDB region_db;
	RegionDB::MapElementIDToRegionID map;
	region_db.read_regions_from_input(i_arr, map);

	EXPECT_EQ(0, ( *(map.find(0)) ).second );
	EXPECT_EQ(1, ( *(map.find(1)) ).second );
	EXPECT_EQ(1, ( *(map.find(3)) ).second );
	EXPECT_EQ(0, ( *(map.find(5)) ).second );
	EXPECT_EQ(1, ( *(map.find(8)) ).second );
	EXPECT_EQ(0, ( *(map.find(9)) ).second );
}


//**************************************************************************************************

void init_db(RegionDB &db,int bc_size, int bulk_size) {
    int i;
    for(i=0; i<bc_size; i++) {
        db.add_region(i, boost::lexical_cast<std::string>(i),1);
    }
    for(; i<bc_size+bulk_size; i++) {
        db.add_region(i, "."+boost::lexical_cast<std::string>(i),1);
    }
}

struct Item {
    unsigned int id;
    std::string l;
    unsigned int dim;
};

void init_map(std::map<unsigned int, Item> &map,unsigned int size) {

    for(unsigned int i=0; i<size; i++) {
        Item xx={i, "xyz", 1};
        map[i]=xx;
    }
}

/***
 * Speed measurments
 *
 * 1M steps
 * Debug    10      DB      1040
 * Debug    10      map     250
 * Debug    100     DB      991
 * Debug    100     map     322
 *
 * 10M steps
 * O3       10      DB      331
 * O3       10      map     132
 * O3       100     DB      330
 * O3       100     map     164
 * 10M steps
 * O3       10      DBmap      452
 * O3       100     DBmap      580
 *
 * O3       100     add_region   4153
 * O3       100     add_region_consistancy_check  825
 * O3       100     add_region_consistancy_check && using iterators  446
 */
/*
#define STEPS (10*1000*1000)

// RegionDB add_item(id, dim) overhead.
TEST(RegionDB, speed_add_region_id) {
        RegionDB region_db;
        init_db(region_db,50,50);
        int ii=0;

        for(int step=0;step < STEPS; step++) {
            ii+= region_db.add_region(3, 1).idx();
            ii+= region_db.add_region(9, 1).idx();
        }
   cout << ii << endl;
}

// boost multi index overhead
TEST(RegionDB, speed_find_id) {
        RegionDB region_db;
        init_db(region_db,50,50);
        int ii=0;

        for(int step=0;step < STEPS; step++) {
            ii+= region_db.find_id(3).idx();
            ii+= region_db.find_id(9).idx();
        }
   cout << ii << endl;
}

// Simplest implementation for ID lookup.
TEST(RegionDB, speed_map) {
        std::map<unsigned int, Item> map;
        init_map(map, 100);
        int ii=0;

        for(int step=0;step < STEPS; step++) {
            ii+= map.find(3)->second.id;
            ii+= map.find(9)->second.id;
        }
   cout << ii << endl;
}
*/
