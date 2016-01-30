/*
 * material_dispatch_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "mesh/region.hh"
#include "mesh/mesh.h"
#include "mesh/region_set.hh"
#include "input/type_base.hh"
#include "input/type_output.hh"
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"
#include <map>

#include <boost/lexical_cast.hpp>

TEST(Region, all) {
	::testing::FLAGS_gtest_death_test_style = "threadsafe";

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
	Input::Type::Array region_set_array_input_type( RegionDB::get_region_set_input_type() );
	Input::ReaderToStorage json_reader( read_sets_json, region_set_array_input_type, Input::FileFormat::format_JSON);
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

	Input::Type::Array element_map_array_input_type( RegionDB::get_region_input_type() );
	Input::ReaderToStorage json_reader( read_element_map_json, element_map_array_input_type, Input::FileFormat::format_JSON);
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


const string read_new_regions = R"JSON(
[
	{
		TYPE="From_Id",
		name = "label_0",
		id = 0
	},
	{
		TYPE="From_Elements",
		name = "label_1",
		id = 1,
		element_list = [8, 5, 6]
	},
	{
		TYPE="From_Label",
		name = "label_2",
		mesh_label = "3D back"
	},
	{
		TYPE="Union",
		name = "label_3",
		regions = ["label_0", "label_1"]
	},
	{
		TYPE="Union",
		name = "label_4",
		regions = ["label_0", "label_2"]
	},
	{
		TYPE="Union",
		name = "label_5",
		region_ids = [37, 38],
		regions = ["label_0", "label_1", "label_2"]
	},
	{
		TYPE="Difference",
		name = "label_6",
		regions = ["label_3", "label_4"]
	},
	{
		TYPE="Intersection",
		name = "label_7",
		regions = ["label_3", "label_4", "label_5"]
	}
]
)JSON";

TEST(NewRegion, read_new_regions) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Profiler::initialize();

    FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
    Mesh * mesh = new Mesh;
    ifstream in(string( mesh_file ).c_str());
    mesh->read_gmsh_from_stream(in, false);

	Input::Type::Array element_map_array_input_type( RegionSetBase::get_input_type() );
	Input::ReaderToStorage json_reader( read_new_regions, element_map_array_input_type, Input::FileFormat::format_JSON);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	mesh->read_regions_from_input(i_arr);

	const RegionDB & region_db = mesh->region_db();
	region_db.print_region_table(cout);

	EXPECT_EQ( 1, region_db.get_region_set("label_0").size() );
	EXPECT_EQ( 0, region_db.get_region_set("label_0")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_1").size() );
	EXPECT_EQ( 1, region_db.get_region_set("label_1")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_2").size() );
	EXPECT_EQ(39, region_db.get_region_set("label_2")[0].id() );

	EXPECT_EQ( 2, region_db.get_region_set("label_3").size() );
	EXPECT_EQ( 0, region_db.get_region_set("label_3")[0].id() );
	EXPECT_EQ( 1, region_db.get_region_set("label_3")[1].id() );

	EXPECT_EQ( 2, region_db.get_region_set("label_4").size() );
	EXPECT_EQ(39, region_db.get_region_set("label_4")[0].id() );
	EXPECT_EQ( 0, region_db.get_region_set("label_4")[1].id() );

	EXPECT_EQ( 5, region_db.get_region_set("label_5").size() );
	EXPECT_EQ(37, region_db.get_region_set("label_5")[0].id() );
	EXPECT_EQ(38, region_db.get_region_set("label_5")[1].id() );
	EXPECT_EQ(39, region_db.get_region_set("label_5")[2].id() );
	EXPECT_EQ( 0, region_db.get_region_set("label_5")[3].id() );
	EXPECT_EQ( 1, region_db.get_region_set("label_5")[4].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_6").size() );
	EXPECT_EQ( 1, region_db.get_region_set("label_6")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_7").size() );
	EXPECT_EQ( 0, region_db.get_region_set("label_7")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("3D back").size() );
	EXPECT_EQ(39, region_db.get_region_set("3D back")[0].id() );

	EXPECT_EQ(37, mesh->element[0].region().id() );
	EXPECT_EQ(39, mesh->element[3].region().id() );
	EXPECT_EQ( 1, mesh->element[4].region().id() );
	EXPECT_EQ( 1, mesh->element[5].region().id() );
	EXPECT_EQ( 1, mesh->element[7].region().id() );
	EXPECT_EQ("label_2", mesh->element[3].region().label() );
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
#ifdef FLOW123D_RUN_UNIT_BENCHMARKS
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
#endif // FLOW123D_RUN_UNIT_BENCHMARKS
