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


class BIHTree_test : public BIHTree {
public:
	BIHTree_test(Mesh* mesh, unsigned int soft_leaf_size_limit)
	: BIHTree(mesh, soft_leaf_size_limit) {}

	/// Tests basic tree parameters (depths, counts of elements)
	void test_tree_params() {
		unsigned int sum_depth = 0;
		unsigned int max_depth = 0;
		unsigned int min_depth = 32767;
		unsigned int leaf_nodes = 0;
		vector<unsigned int> elements; // counts of elements in leaf nodes

		for (unsigned int i=0; i<nodes_.size(); i++) {
			if (nodes_[i].is_leaf()) {
				if (nodes_[i].depth() > max_depth) max_depth = nodes_[i].depth();
				if (nodes_[i].depth() < min_depth) min_depth = nodes_[i].depth();
				sum_depth += nodes_[i].depth();
				elements.push_back(nodes_[i].leaf_size());
				++leaf_nodes;
			}
		}

		EXPECT_EQ(elements.size(), leaf_nodes);
		std::sort(elements.begin(), elements.end());

		double avg_depth = (double) sum_depth / (double) leaf_nodes;
		double median_elements = (double)(elements[ leaf_nodes/2 ] + elements[ (leaf_nodes-1)/2 ]) / 2.0;

		cout << endl << "-------------------------";
		cout << endl << "BIH tree parameters:" << endl;
		cout << "- maximal depth: " << max_depth;
		cout << ", minimal depth: " << min_depth;
		cout << ", average depth: " << avg_depth << endl;
		cout << "- median of elements in leaf nodes: " << median_elements << endl;
		cout << "- minimal elements in leaf nodes: " << elements[ 0 ] << endl;
		cout << "- maximal elements in leaf nodes: " << elements[ elements.size()-1 ] << endl;
	}

	/// Printout structure of BIH tree
	void BIH_output() {
		cout << endl << "-------------------------";
		cout << endl << "BIH tree output:";
		BIH_output_node(0);
		cout << endl;
	}

protected:
	/**
	 * Returns count of elements of get node
	 *  - count of elements in leaf node
	 *  - sum of elements of all leaf descendants for inner node
	 */
	unsigned int BIH_elements_in_node(unsigned int node_index) {
		const BIHNode &node = nodes_[node_index];
		if (node.is_leaf()) {
			return node.leaf_size();
		} else {
			return BIH_elements_in_node( node.child(0) ) + BIH_elements_in_node( node.child(1) );
		}
	}

	/// Printout node of BIH tree
	void BIH_output_node(unsigned int node_index, unsigned int depth = 0) {
		const BIHNode &node = nodes_[node_index];

		cout << endl;
		if (depth>1) for (unsigned int i=0; i<depth-1; i++) cout << "| ";
		if (node.is_leaf()) {
			cout << "- leaf node: idx = " << node_index << ", " << node.leaf_size() << " elements, depth " << (unsigned int)node.depth();
		} else {
			if (depth>0) cout << "| ";
			cout << "inner node: idx = " << node_index << ", " << BIH_elements_in_node(node_index) << " elements, children idx (";
			cout << node.child(0) << "," << node.child(1) << ")";
			BIH_output_node( node.child(0), depth+1 );
			BIH_output_node( node.child(1), depth+1 );
		}
	}
};


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
	    bt = new BIHTree_test(mesh, leaf_size_limit);
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

		bt->test_tree_params();
		//bt->BIH_output();
		cout << endl;
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
			for(unsigned int j=0; j< bf_result.size(); j++) {
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
			for(unsigned int j=0; j< bf_point_result.size(); j++) {
				EXPECT_EQ(bf_point_result[j], result_point_vec[j]);
			}
		}
	}


	BIH_test()
	: r_gen(123), mesh(nullptr), bt(nullptr)
	{
        Profiler::initialize();
	}

	~BIH_test() {
		if (mesh !=nullptr) delete mesh;
		if (bt !=nullptr) delete bt;
		mesh = nullptr;
		bt = nullptr;

		Profiler::uninitialize();
	}

	std::mt19937	r_gen;
	Mesh *mesh;
	BIHTree_test *bt;
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

/**
 * Unit test of BIH tree on large mesh (111 000 elements).
 *
 * Only geo file is committed to repository. Mesh file is stored
 * on 'http://bacula.nti.tul.cz/~jan.brezina/BIH_refined.tar.gz'
 * in directory '1024_111324_el'
 */
/*TEST_F(BIH_test, find_bounding_box_6) {
	FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_111324_elem.msh", FilePath::input_file);
	this->create_tree(mesh_file);

	this->test_find_boxes();
}*/


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

