/*
 * mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#include <gtest/gtest.h>
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include <iostream>
#include <vector>

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

    Mesh mesh;
    GmshMeshReader reader(mesh_file);
    reader.read_mesh(&mesh);

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
