/*
 * mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#define TEST_USE_MPI
#include <gtest_mpi.hh>

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include <iostream>
#include <vector>
#include "mesh/accessors.hh"
#include "input/json_to_storage.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"



using namespace std;

class MeshTest :  public testing::Test, public Mesh {
    virtual void SetUp() {}
    virtual void TearDown() {}
};


TEST_F(MeshTest, intersect_nodes_lists) {
    node_elements.resize(3);
    node_elements[0]={ 0, 1, 2, 3, 4};
    node_elements[1]={ 0, 2, 3, 4};
    node_elements[2]={ 0, 1, 2, 4};

    vector<unsigned int> node_list={0,1,2};
    vector<unsigned int> result;
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,2,4} ), result );

    node_list={0,1};
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,2,3,4} ), result );

    node_list={0};
    intersect_element_lists(node_list, result);
    EXPECT_EQ( vector<unsigned int>( {0,1,2,3,4} ), result );

}



TEST(MeshTopology, make_neighbours_and_edges) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);

    Profiler::initialize();
    
    Mesh mesh;
    ifstream in(string(mesh_file).c_str());
    mesh.read_gmsh_from_stream(in);


    EXPECT_EQ(9, mesh.n_elements());
    EXPECT_EQ(18, mesh.bc_elements.size());

    // check bc_elements
    EXPECT_EQ(101 , mesh.bc_elements[0].region().id() );
    EXPECT_EQ(101 , mesh.bc_elements[1].region().id() );
    EXPECT_EQ(102 , mesh.bc_elements[2].region().id() );
    EXPECT_EQ(102 , mesh.bc_elements[3].region().id() );
    EXPECT_EQ( -3 , int( mesh.bc_elements[4].region().id() ) );
    EXPECT_EQ( -3 , int( mesh.bc_elements[17].region().id() ) );

    //check edges
    EXPECT_EQ(28,mesh.n_edges());

    //check neighbours
    EXPECT_EQ(6, mesh.n_vb_neighbours() );

}


// Test input for mesh
const string mesh_input = R"JSON(
{ 
  mesh_file="mesh/simplest_cube.msh",
  regions=[{id=3000, name="new region", element_list=[6,7]},
           {id=37, name="1D rename"}],
  sets=[{name="3D", region_ids=[39,40]}]
}
)JSON";

TEST(Mesh, init_from_input) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::JSONToStorage reader;
    std::stringstream in(mesh_input);
    reader.read_stream( in, Mesh::input_type );
    Mesh mesh( reader.get_root_interface<Input::Record>() );
    mesh.init_from_input();

    EXPECT_EQ( 37, mesh.element_accessor(0).region().id() );
    EXPECT_EQ( "1D rename", mesh.element_accessor(0).region().label() );

    EXPECT_EQ( 38, mesh.element_accessor(1).region().id() );
    EXPECT_EQ( 38, mesh.element_accessor(2).region().id() );
    EXPECT_EQ( 39, mesh.element_accessor(3).region().id() );
    EXPECT_EQ( 39, mesh.element_accessor(4).region().id() );
    EXPECT_EQ( 3000, mesh.element_accessor(5).region().id() );
    EXPECT_EQ( 3000, mesh.element_accessor(6).region().id() );
    EXPECT_EQ( 40, mesh.element_accessor(7).region().id() );
    EXPECT_EQ( 40, mesh.element_accessor(8).region().id() );

    RegionSet set = mesh.region_db().get_region_set("3D");
    EXPECT_EQ( 39, set[0].id() );
    EXPECT_EQ( 40, set[1].id() );

}
