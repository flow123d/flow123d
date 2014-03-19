/*
 * bih_tree_test.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: jb
 */



#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
#include <cmath>
#include <algorithm>
#include <fstream>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"


class BIH_test : public testing::Test {
public:
	void create_tree(const FilePath &mesh_file) {

		START_TIMER("create mesh");
		GmshMeshReader reader(mesh_file);
		mesh = new Mesh();
		reader.read_mesh(mesh);
		END_TIMER("create mesh");

	    int leaf_size_limit = 10;
	    START_TIMER("create bih tree");
	    bt = new BIHTree(mesh, leaf_size_limit);
	    END_TIMER("create bih tree");

	    EXPECT_EQ(mesh->n_elements(), bt->get_element_count());

	    this->init_random_gen();

	}

	void init_random_gen() {
		auto tbox_min = bt->tree_box().min();
		auto tbox_max = bt->tree_box().max();
		rx_ = std::uniform_real_distribution<>( tbox_min[0] , tbox_max[0] );
		ry_ = std::uniform_real_distribution<>( tbox_min[0] , tbox_max[0] );
		rz_ = std::uniform_real_distribution<>( tbox_min[0] , tbox_max[0] );
	}

	BoundingBox::Point r_point() {
		return BoundingBox::Point( {rx_(this->r_gen), ry_(this->r_gen), rz_(this->r_gen)} );
	}

	void test_find_boxes() {
		test_insec_elements();
		test_insec_points();

		Profiler::instance()->output(MPI_COMM_WORLD, cout);
	}


	void test_insec_elements() {
		for(int i=0; i < n_test_trials; i++) {
			BoundingBox box( vector<BoundingBox::Point>({r_point(), r_point()}) );

			//cout << "\n-------------------------" << endl;
			//cout << "box: " << box <<endl;

			vector<unsigned int> bf_result;
			FOR_ELEMENTS(mesh, ele) {
				EXPECT_EQ( box.intersect(ele->bounding_box()) , ele->bounding_box().intersect(box) );
				if (box.intersect(ele->bounding_box()) ) bf_result.push_back(ele.index());
			}

			vector<unsigned int> result_vec;
			START_TIMER("find bounding box");
			bt->find_bounding_box(box, result_vec);
			END_TIMER("find bounding box");
			std::sort(result_vec.begin(), result_vec.end());

			//cout << endl << "full search: " << endl;
			//for(unsigned int i_el : bf_result) cout << " " << this->mesh->element(i_el).id();
			//cout << endl << "bih search: " << endl;
			//for(unsigned int i_el : result_vec) cout << " " << this->mesh->element(i_el).id();

			ASSERT_EQ(bf_result.size(), result_vec.size());
			for(int j=0; j< bf_result.size(); j++) {
				EXPECT_EQ(bf_result[j], result_vec[j]);
			}
		}
	}


	void test_insec_points() {
		for(int i=0; i < n_test_trials; i++) {
			BoundingBox::Point point( r_point() );

			//cout << "\n-------------------------" << endl;
			//cout << "point: " << point << endl;

			vector<unsigned int> bf_point_result;
			FOR_ELEMENTS(mesh, ele) {
				if (ele->bounding_box().contains_point(point) ) bf_point_result.push_back(ele.index());
			}

			vector<unsigned int> result_point_vec;
			START_TIMER("find point");
			bt->find_point(point, result_point_vec);
			END_TIMER("find point");
			std::sort(result_point_vec.begin(), result_point_vec.end());

			//cout << endl << "full search: " << endl;
			//for(unsigned int i_el : bf_point_result) cout << " " << this->mesh->element(i_el).id();
			//cout << endl << "bih search: " << endl;
			//for(unsigned int i_el : result_point_vec) cout << " " << this->mesh->element(i_el).id();

			ASSERT_EQ(bf_point_result.size(), result_point_vec.size());
			for(int j=0; j< bf_point_result.size(); j++) {
				EXPECT_EQ(bf_point_result[j], result_point_vec[j]);
			}
		}
	}


	BIH_test()
	: mesh(nullptr), bt(nullptr), r_gen(123)
	{
        Profiler::initialize();
	}

	~BIH_test() {
		if (mesh !=nullptr) delete mesh;
		if (bt !=nullptr) delete bt;
		mesh = nullptr;
		bt = nullptr;
	}

	std::mt19937	r_gen;
	Mesh *mesh;
	BIHTree *bt;
	std::uniform_real_distribution<> rx_, ry_, rz_;
	const static int n_test_trials=5;
};


TEST_F(BIH_test, find_bounding_box_1) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_108_elem.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}

TEST_F(BIH_test, find_bounding_box_2) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/fields/simplest_cube_3d.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}

TEST_F(BIH_test, find_bounding_box_3) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_7590_elem.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}

TEST_F(BIH_test, find_bounding_box_4) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_188_elem.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}

TEST_F(BIH_test, find_bounding_box_5) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_27936_elem.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}

/*
/// Generates random double number in interval <fMin, fMax>
double f_rand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


/// Gets count of intersected elements with bounding box
unsigned int get_intersection_count(BoundingBox &bb, std::vector<BoundingBox> &boundingBoxes) {
	unsigned int insecElements = 0;

	for (unsigned int i=0; i<boundingBoxes.size(); i++) {
		if (bb.intersect(boundingBoxes[i])) insecElements++;
	}

	return insecElements;
}

*/
/**
 * Creates tree and performs its tests
 *  - creates tree from mesh file
 *  - performs tests of basic parameters (maximal depth, count of nodes)
 *  - tests intersection with bounding box out of mesh
 *  - tests intersection with three bounding boxes in mesh
 */
/*
void create_test_tree(FilePath &meshFile, unsigned int elementLimit = 20) {
        Profiler::initialize();
	unsigned int maxDepth, minDepth, sumDepth, leafNodesCount, innerNodesCount, sumElements, insecSize;
	double avgDepth;
	Mesh mesh;
	GmshMeshReader reader(meshFile);
	BoundingBox bb;
	arma::vec3 min, max;
	std::vector<unsigned int> searchedElements;

	reader.read_mesh(&mesh);

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
	bb=BoundingBox(arma::vec3("0 0 1.01"), arma::vec3("0.1 0.1 1.05"));
	bt.find_bounding_box(bb, searchedElements);
	EXPECT_EQ(0, searchedElements.size());

	// tests of intersection with bounding box in mesh near point [-1, -1, -1]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(-0.99, -0.97);
		max(i) = f_rand(-0.96, -0.94);
	}
	bb=BoundingBox(min, max);
	bt.find_bounding_box(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements()); // get intersections by linear search
	EXPECT_EQ(searchedElements.size(), insecSize);

	// tests of intersection with bounding box in mesh near point [0, 0, 0]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(-0.03, -0.01);
		max(i) = f_rand(+0.01, +0.03);
	}
	bb=BoundingBox(min, max);
	bt.find_bounding_box(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements());
	EXPECT_EQ(searchedElements.size(), insecSize);

	// tests of intersection with bounding box in mesh near point [0.1, 0.5, 0.9]
	for (int i=0; i<3; i++) {
		min(i) = f_rand(0.07 + i * 0.4, 0.09 + i * 0.4);
		max(i) = f_rand(0.11 + i * 0.4, 0.13 + i * 0.4);
	}
	bb=BoundingBox(min, max);
	bt.find_bounding_box(bb, searchedElements);
	insecSize = get_intersection_count(bb, bt.get_elements());
	EXPECT_EQ(searchedElements.size(), insecSize);

	// tests of intersection with point in mesh near point [0.2, 0.3, 0.4]
	insecSize = 0;
	for (int i=0; i<3; i++) {
		min(i) = f_rand(0.18 + i * 0.1, 0.22 + i * 0.1);
	}
	bt.find_point(min, searchedElements);
	for (unsigned int i=0; i<bt.get_elements().size(); i++) {
		if ( bt.get_elements()[i].contains_point(min) ) insecSize++;
	}
	EXPECT_EQ(searchedElements.size(), insecSize);
}

/*
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



*/


TEST(BIH_Tree_Test, 2d_mesh) {
    Profiler::initialize();
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/noncompatible_small.msh", FilePath::input_file);

    Mesh mesh;
	GmshMeshReader reader(mesh_file);

	reader.read_mesh(&mesh);
	unsigned int element_limit=20;
	BIHTree bt(&mesh, element_limit);
	std::vector<unsigned int> insec_list;

	bt.find_bounding_box(BoundingBox(arma::vec3("-1.1 0 0"), arma::vec3("-0.7 0 0")), insec_list);
	for(auto i_ele : insec_list) {
		cout << "idx: " << i_ele << "id: " << mesh.element.get_id( &(mesh.element[i_ele]) ) << endl;
	}


}

