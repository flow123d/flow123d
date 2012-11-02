/*
 * bih_tree_test.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: jb
 */




/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define DEBUG

#include <gtest/gtest.h>
#include <sstream>
#include <string>

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bih_tree.hh"

// simplest cube 123d
string gmsh_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
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
9
1 1 2 37 20 7 3
2 2 2 38 34 6 3 7
3 2 2 38 36 3 1 7
4 4 2 39 40 3 7 1 2
5 4 2 39 40 3 7 2 8
6 4 2 39 40 3 7 8 6
7 4 2 39 42 3 7 6 5
8 4 2 39 42 3 7 5 4
9 4 2 39 42 3 7 4 1
$EndElements
)CODE";


void create_tree(FilePath &meshFile, int elementLimit = 0) {
	int maxDepth, minDepth, sum, leaves;
	Mesh mesh;
	GmshMeshReader reader;

	reader.read(meshFile, &mesh);

	BIHTree bt(&mesh, elementLimit);
	bt.get_tree_depth(maxDepth, minDepth, sum, leaves, false);
}

TEST(BIHTree_Test, mesh_from_stream) {
    stringstream ss(gmsh_mesh);

    Mesh mesh;
    GmshMeshReader reader;

    reader.read(ss, &mesh);

    BIHTree bt(&mesh, 4);
    int maxDepth, minDepth, sum, leaves;
    bt.get_tree_depth(maxDepth, minDepth, sum, leaves, false);
}

TEST(BIHTree_Test, mesh_108_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_108_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_390_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_390_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_1907_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_1907_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_7590_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_7590_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_31949_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_31949_elem.msh", FilePath::input_file);

    create_tree(mesh_file, 100);
}

TEST(BIHTree_Test, mesh_188_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_188_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_482_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_482_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_1638_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_1638_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_5927_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_5927_elem.msh", FilePath::input_file);

    create_tree(mesh_file);
}

TEST(BIHTree_Test, mesh_27936_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_27936_elem.msh", FilePath::input_file);

    create_tree(mesh_file, 100);
}

TEST(BIHTree_Test, mesh_111324_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_111324_elem.msh", FilePath::input_file);

    create_tree(mesh_file, 100);
}


