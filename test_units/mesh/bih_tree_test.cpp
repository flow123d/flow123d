/*
 * bih_tree_test.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: jb
 */


#define DEBUG

#include <gtest/gtest.h>
#include <cmath>

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bih_tree.hh"


/// Generates random double number in interval <fMin, fMax>
double f_rand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


/// Gets count of intersected elements with bounding box
int get_intersection_count(BoundingBox &bb, std::vector<BoundingBox> &boundingBoxes) {
	int insecElements = 0;

	for (int i=0; i<boundingBoxes.size(); i++) {
		if (bb.intersection(boundingBoxes[i])) insecElements++;
	}

	return insecElements;
}


/**
 * Creates tree and performs its tests
 *  - creates tree from mesh file
 *  - performs tests of basic parameters (maximal depth, count of nodes)
 *  - tests intersection with bounding box out of mesh
 *  - tests intersection with three bounding boxes in mesh
 */
void create_test_tree(FilePath &meshFile, int elementLimit = 20) {
	int maxDepth, minDepth, sumDepth, leafNodesCount, innerNodesCount, sumElements, insecSize;
	double avgDepth;
	Mesh mesh;
	GmshMeshReader reader;
	BoundingBox bb;
	arma::vec3 min, max;
	std::vector<int> searchedElements;

	reader.read(meshFile, &mesh);

	// creates tree and tests its basic parameters
	BIHTree bt(&mesh, elementLimit);
	bt.get_tree_params(maxDepth, minDepth, avgDepth, leafNodesCount, innerNodesCount, sumElements);
	// For ideal case the node count should be comparable to total number of elements devided
	// by number of elements in leaf nodes. In this test we check that the overhead is not
	// too big. Empirical test showed, that for big meshes the actual number of nodes is
	// 3.7 * ideal number of nodes
	EXPECT_LT( (leafNodesCount + innerNodesCount) , 4 * (mesh.n_elements() / elementLimit) );
	// Check  the conditions that limit growth of the redundancy of elements in nodes
	// The condition says that number of elements in the child node has to have number of elements less then
	// 0.8 * (number of elements in parent node).
	EXPECT_LT( maxDepth, (log2(elementLimit) - log2(mesh.n_elements())) / log2(0.8) );

	// tests of intersection with bounding box out of mesh
	bb.set_bounds(arma::vec3("0 0 1.01"), arma::vec3("0.1 0.1 1.05"));
	bt.find_elements(bb, searchedElements);
	EXPECT_EQ(0, searchedElements.size());

	// tests of intersection with bounding box in mesh near point [-1, -1, -1]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(-0.99, -0.97);
		max(i) = f_rand(-0.96, -0.94);
	}
	bb.set_bounds(min, max);
	bt.find_elements(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements()); // get intersections by linear search
	EXPECT_EQ(searchedElements.size(), insecSize);

	// tests of intersection with bounding box in mesh near point [0, 0, 0]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(-0.03, -0.01);
		max(i) = f_rand(+0.01, +0.03);
	}
	bb.set_bounds(min, max);
	bt.find_elements(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements());
	EXPECT_EQ(searchedElements.size(), insecSize);

	// tests of intersection with bounding box in mesh near point [0.1, 0.5, 0.9]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(0.07 + i * 0.4, 0.09 + i * 0.4);
		max(i) = f_rand(0.11 + i * 0.4, 0.13 + i * 0.4);
	}
	bb.set_bounds(min, max);
	bt.find_elements(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements());
	EXPECT_EQ(searchedElements.size(), insecSize);

}


TEST(BIHTree_Test, mesh_108_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_108_elem.msh", FilePath::input_file);

    create_test_tree(mesh_file);

    // tree->number_of_nodes  < C * (mesh.n_elements() / areaElementLimit) -- C=4
    // tree->max_depth  < (lg(A) - lg(N))/lg(0.8)
    //
    // find_elements - box mimo sit - EXPECT_EQ( 0, result.size())
    // find_elements - box uvnitr - porovnat s uplnym pruchodem site
    // asi tri boxy .. EXPECT_EQ

}


TEST(BIHTree_Test, mesh_7590_elements_homogeneous) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_7590_elem.msh", FilePath::input_file);

    create_test_tree(mesh_file);
}


TEST(BIHTree_Test, mesh_188_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_188_elem.msh", FilePath::input_file);

    create_test_tree(mesh_file);
}


TEST(BIHTree_Test, mesh_27936_elements_refined) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_27936_elem.msh", FilePath::input_file);

    create_test_tree(mesh_file);
}

