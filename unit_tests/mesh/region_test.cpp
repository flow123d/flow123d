/*
 * material_dispatch_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest.hh>

#include "mesh/region.hh"
#include "mesh/mesh.h"
#include "mesh/region_set.hh"
#include "mesh/msh_gmshreader.h"
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
    Region r=region_db.add_region(1003, region_db.create_label_from_id(1003), 3);
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
		Region r=region_db.add_region(1, "user_defined_name", RegionDB::undefined_dim);
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


// simplest mesh
string gmsh_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
6
1       37      "1D diagonal"
2       38      "2D XY diagonal"
2       101     ".top side"
2       102     ".bottom side"
3       39      "3D back"
3       40      "3D front"
$EndPhysicalNames
$Nodes
8
1 1 1 1
2 -1 1 1
3 -1 -1 1
4 1 -1 1
5 1 -1 -1
6 -1 -1 -1
7 1 1 -1
8 -1 1 -1
$EndNodes
$Elements
13
1 1 2 37 20 7 3
2 2 2 38 34 6 3 7
3 2 2 38 36 3 1 7
4 4 2 39 40 3 7 1 2
5 4 2 39 40 3 7 2 8
6 4 2 39 40 3 7 8 6
7 4 2 40 42 3 7 6 5
8 4 2 40 42 3 7 5 4
9 4 2 41 42 3 7 4 1
10 2 2 101 101 1 2 3
11 2 2 101 101 1 3 4
12 2 2 102 102 6 7 8
13 2 2 102 102 7 6 5 
$EndElements
)CODE";

const string read_regions_yaml = R"YAML(
- !From_Id
  name: 3D front rename
  id: 40
- !From_Label
  name: 3D back rename
  mesh_label: 3D back
- !From_Elements
  name: label_0
  element_list:
   - 6
- !From_Elements
  name: label_1
  id: 1
  element_list:
   - 8
   - 5
- !Union
  name: label_2
  regions:
   - 3D front rename
   - label_1
- !Union
  name: label_3
  regions: 
   - 3D front rename
   - 3D back rename
- !Union
  name: label_4
  region_ids: 
   - 37
   - 38
   - 41
  regions:
   - 3D front rename
   - label_1
   - 3D back rename
- !Difference
  name: label_5
  regions:
   - label_2
   - label_3
- !Intersection
  name: label_6
  regions:
   - label_2
   - label_3
   - label_4
)YAML";

TEST(Region, read_regions_from_yaml) {
    Profiler::initialize();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    stringstream in(gmsh_mesh.c_str());
    GmshMeshReader reader(in);
    Mesh * mesh = new Mesh;

	Input::Type::Array element_map_array_input_type( RegionSetBase::get_input_type() );
	Input::ReaderToStorage json_reader( read_regions_yaml, element_map_array_input_type, Input::FileFormat::format_YAML);
	Input::Array i_arr = json_reader.get_root_interface<Input::Array>();

	reader.read_physical_names(mesh);
	mesh->read_regions_from_input(i_arr);
	reader.read_mesh(mesh);
	mesh->check_and_finish();

	const RegionDB & region_db = mesh->region_db();
	region_db.print_region_table(cout);

	EXPECT_EQ( 1, region_db.get_region_set("3D front rename").size() );
	EXPECT_EQ(40, region_db.get_region_set("3D front rename")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("3D back rename").size() );
	EXPECT_EQ(39, region_db.get_region_set("3D back rename")[0].id() );

	EXPECT_EQ(  1, region_db.get_region_set("label_0").size() );
	EXPECT_EQ(103, region_db.get_region_set("label_0")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_1").size() );
	EXPECT_EQ( 1, region_db.get_region_set("label_1")[0].id() );

	EXPECT_EQ( 2, region_db.get_region_set("label_2").size() );
	EXPECT_EQ(40, region_db.get_region_set("label_2")[0].id() );
	EXPECT_EQ( 1, region_db.get_region_set("label_2")[1].id() );

	EXPECT_EQ( 2, region_db.get_region_set("label_3").size() );
	EXPECT_EQ(39, region_db.get_region_set("label_3")[0].id() );
	EXPECT_EQ(40, region_db.get_region_set("label_3")[1].id() );

	EXPECT_EQ( 6, region_db.get_region_set("label_4").size() );
	EXPECT_EQ(37, region_db.get_region_set("label_4")[0].id() );
	EXPECT_EQ(38, region_db.get_region_set("label_4")[1].id() );
	EXPECT_EQ(39, region_db.get_region_set("label_4")[2].id() );
	EXPECT_EQ(40, region_db.get_region_set("label_4")[3].id() );
	EXPECT_EQ( 1, region_db.get_region_set("label_4")[4].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_5").size() );
	EXPECT_EQ( 1, region_db.get_region_set("label_5")[0].id() );

	EXPECT_EQ( 1, region_db.get_region_set("label_6").size() );
	EXPECT_EQ(40, region_db.get_region_set("label_6")[0].id() );

	EXPECT_EQ( 9, region_db.get_region_set("ALL").size() );
	EXPECT_EQ( 7, region_db.get_region_set("BULK").size() );
	EXPECT_EQ( 2, region_db.get_region_set(".BOUNDARY").size() );

	EXPECT_EQ( 37, mesh->element[0].region().id() );
	EXPECT_EQ( 39, mesh->element[3].region().id() );
	EXPECT_EQ(  1, mesh->element[4].region().id() );
	EXPECT_EQ(103, mesh->element[5].region().id() );
	EXPECT_EQ( 40, mesh->element[6].region().id() );
	EXPECT_EQ(  1, mesh->element[7].region().id() );
	EXPECT_EQ("label_1", mesh->element[4].region().label() );
}


const string invalid_input_1 = R"YAML(
- !From_Id
  name: 3D front
  id: 37
)YAML";

const string invalid_input_2 = R"YAML(
- !From_Elements
  name: 3D back
  id: 41
  element_list:
   - 6
)YAML";

const string invalid_input_3 = R"YAML(
- !Intersection
  name: insec
  regions:
   - 3D front
   - 3D back
)YAML";

TEST(Region, read_regions_error_messages) {
    Profiler::initialize();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
    GmshMeshReader reader(mesh_file);
    Mesh * mesh = new Mesh;
	reader.read_physical_names(mesh);

	{
		Input::Type::Array element_map_array_input_type( RegionSetBase::get_input_type() );
		Input::ReaderToStorage json_reader( invalid_input_1, element_map_array_input_type, Input::FileFormat::format_YAML);
		Input::Array i_arr = json_reader.get_root_interface<Input::Array>();
		EXPECT_THROW_WHAT( { mesh->read_regions_from_input(i_arr); }, RegionDB::ExcNonuniqueLabel, "id: 37, label: '3D front'");
	}
	{
		Input::Type::Array element_map_array_input_type( RegionSetBase::get_input_type() );
		Input::ReaderToStorage json_reader( invalid_input_2, element_map_array_input_type, Input::FileFormat::format_YAML);
		Input::Array i_arr = json_reader.get_root_interface<Input::Array>();
		EXPECT_THROW_WHAT( { mesh->read_regions_from_input(i_arr); }, RegionDB::ExcNonuniqueLabel, "id: 41, label: '3D back'");
	}
	{
		Input::Type::Array element_map_array_input_type( RegionSetBase::get_input_type() );
		Input::ReaderToStorage json_reader( invalid_input_3, element_map_array_input_type, Input::FileFormat::format_YAML);
		Input::Array i_arr = json_reader.get_root_interface<Input::Array>();
		EXPECT_THROW_WHAT( { mesh->read_regions_from_input(i_arr); }, RegionSetBase::ExcEmptyRegionSetResult,
				"Empty result of Intersection operation.");
	}
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
TEST(RegionDB, speed_get_region_id) {
        RegionDB region_db;
        init_db(region_db,50,50);
        int ii=0;

        for(int step=0;step < STEPS; step++) {
            ii+= region_db.get_region(3, 1).idx();
            ii+= region_db.get_region(9, 1).idx();
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
