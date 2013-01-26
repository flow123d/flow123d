/*
 * mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#include <gtest/gtest.h>
#include "mesh/mesh.h"
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
