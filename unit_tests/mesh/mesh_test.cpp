/*
 * mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include <iostream>
#include <vector>
#include "mesh/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"



using namespace std;

class MeshTest :  public testing::Test, public Mesh {
public:
    MeshTest()
    : Mesh(*mesh_constructor())
    {
    }

    ~MeshTest()
    {
    }
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
    
    Mesh * mesh = mesh_constructor();
    ifstream in(string(mesh_file).c_str());
    GmshMeshReader gmsh_reader( in );
    gmsh_reader.read_mesh(mesh);

    EXPECT_EQ(9, mesh->n_elements());
    EXPECT_EQ(18, mesh->bc_elements.size());

    // check bc_elements
    EXPECT_EQ(101 , mesh->bc_elements[0].region().id() );
    EXPECT_EQ(101 , mesh->bc_elements[1].region().id() );
    EXPECT_EQ(102 , mesh->bc_elements[2].region().id() );
    EXPECT_EQ(102 , mesh->bc_elements[3].region().id() );
    EXPECT_EQ( -3 , int( mesh->bc_elements[4].region().id() ) );
    EXPECT_EQ( -3 , int( mesh->bc_elements[17].region().id() ) );

    //check edges
    EXPECT_EQ(28,mesh->n_edges());

    //check neighbours
    EXPECT_EQ(6, mesh->n_vb_neighbours() );

    delete mesh;

}


const string mesh_input = R"YAML(
mesh_file: "mesh/simplest_cube.msh"
regions:
 - !From_Elements
   id: 3000
   name: new region
   element_list:
    - 6
    - 7
 - !From_Id
   id: 37
   name: 1D rename
 - !Union
   name: 3D
   region_ids:
    - 39
    - 40
)YAML";

TEST(Mesh, init_from_input) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_constructor(mesh_input, Input::FileFormat::format_YAML);
    GmshMeshReader gmsh_reader( mesh->mesh_file() );
    gmsh_reader.read_mesh(mesh);

    EXPECT_EQ( 37, mesh->element_accessor(0).region().id() );
    EXPECT_EQ( "1D rename", mesh->element_accessor(0).region().label() );

    EXPECT_EQ( 38, mesh->element_accessor(1).region().id() );
    EXPECT_EQ( 38, mesh->element_accessor(2).region().id() );
    EXPECT_EQ( 39, mesh->element_accessor(3).region().id() );
    EXPECT_EQ( 39, mesh->element_accessor(4).region().id() );
    EXPECT_EQ( 3000, mesh->element_accessor(5).region().id() );
    EXPECT_EQ( 3000, mesh->element_accessor(6).region().id() );
    EXPECT_EQ( 40, mesh->element_accessor(7).region().id() );
    EXPECT_EQ( 40, mesh->element_accessor(8).region().id() );

    RegionSet set = mesh->region_db().get_region_set("3D");
    EXPECT_EQ( 39, set[0].id() );
    EXPECT_EQ( 40, set[1].id() );

    delete mesh;
}


// simplest mesh
string small_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 -1 1 1
2 -1 -1 1
3 1 1 -1
4 -1 1 -1
$EndNodes
$Elements
1
1 4 2 39 40 2 3 1 4
$EndElements
)CODE";

TEST(Mesh, decompose_problem) {
    Mesh * mesh = mesh_constructor();
    stringstream in(small_mesh.c_str());
    GmshMeshReader gmsh_reader( in );
    EXPECT_THROW_WHAT( { gmsh_reader.read_mesh(mesh); }, Partitioning::ExcDecomposeMesh,
    		"greater then number of elements 1. Can not make partitioning of the mesh");

    delete mesh;
}
