/**!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: bih_tree.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "mesh/bih_tree.hh"
#include "mesh/bih_node.hh"
#include "system/global_defs.h"
#include <ctime>

/**
 * Minimum reduction of number of elements in child nodes to allow
 * splitting of a node during tree creation.
 */
const double BIHTree::min_reduce_factor = 0.8;
/**
 * Maximum grow factor of total number of elements in childs compared to number of elements in
 * the parent node during tree creation.
 */
const double BIHTree::max_grow_factor = 1.5;


BIHTree::BIHTree(Mesh* mesh, unsigned int soft_leaf_size_limit) {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	srand((unsigned)time(0));

	mesh_ = mesh;
	nodes_.reserve(2 * mesh_->n_elements() / soft_leaf_size_limit);

	//START_TIMER("BIH Tree");

	element_boxes();
	create_tree(soft_leaf_size_limit);

	free_memory();

	//END_TIMER("BIH Tree");

	xprintf(Msg, " - Tree created\n");
}



BIHTree::~BIHTree() {

}


void BIHTree::create_tree(unsigned int soft_leaf_size_limit) {
	create_root_node();

	// create tree
	while (queue_.size()) {
		BIHNode child[BIHNode::child_count];
		BoundingBox child_box[BIHNode::child_count];

		ASSERT(queue_.front() < nodes_.size(), "Idx %d out of nodes_ of size %d.\n", queue_.front(), nodes_.size());

		BIHNode & actual_node = nodes_[queue_.front()];

		ASSERT(actual_node.is_leaf(), "Not leaf: %d\n",queue_.front() );
		unsigned char depth = actual_node.depth();
		unsigned int actual_node_leaf_size = actual_node.leaf_size();

		if (actual_node_leaf_size <= soft_leaf_size_limit) {
			// mark node as definitive leaf
			put_leaf_elements();
			continue;
		}

		unsigned char splitting_axis = box_queue_.front().longest_axis();
		set_node_median(splitting_axis,actual_node);
		actual_node.axis_ = splitting_axis; // actual_node is not leaf anymore


		//calculate bounding boxes of subareas and create them
		box_queue_.front().split(actual_node.axis(),actual_node.median(),child_box[0], child_box[1]);
		child[0].set_depth(depth+1);
		child[1].set_depth(depth+1);

	    distribute_elements(child[0], child[1]);

		// test count of elements in subareas
		if (child[0].leaf_size() > min_reduce_factor * actual_node_leaf_size ||
		    child[1].leaf_size() > min_reduce_factor * actual_node_leaf_size ||
		    child[0].leaf_size() + child[1].leaf_size() > max_grow_factor * actual_node_leaf_size) {

			list_element_index_next_.erase(list_element_index_next_.begin() + child[0].child_[0],
										   list_element_index_next_.begin() + child[1].child_[1]);
			actual_node.set_depth(depth);
			put_leaf_elements();
		    continue;
		}

		// put nodes into vectors
		test_new_level();
		for (unsigned int i=0; i<BIHNode::child_count; i++) {
			// can not use actual_node, since it can be invalid due to reallocation of nodes_

		    nodes_[queue_.front()].child_[i] = nodes_.size();
			queue_.push_back(nodes_.size());
			box_queue_.push_back(child_box[i]);
			nodes_.push_back(child[i]);
		}

		// remove processed node from queue
		queue_.pop_front();
		box_queue_.pop_front();
	}

}


void BIHTree::create_root_node() {

	BoundingBox main_box;

	BIHNode bih_node(0);
	bih_node.child_[0] = 0;
	bih_node.child_[1] = mesh_->n_elements();
	nodes_.push_back(bih_node);

	// put indexes of all elements to vector (for first level of tree)
	list_element_index_.resize(mesh_->n_elements());
	for (unsigned int i=0; i<mesh_->n_elements(); i++) {
		list_element_index_[i] = i;
	}

	// find minimal and maximal coordination of whole mesh
	Node* node = mesh_->node_vector.begin();
	main_box = BoundingBox(node->point(), node->point());

	FOR_NODES(mesh_, node ) {
		main_box.expand( node->point() );
	}

	// put index of root node and its coordinations to vectors
	queue_.clear();
	queue_.push_back(0);
	box_queue_.push_back(main_box);
}




void BIHTree::set_node_median(unsigned char axis, BIHNode &node)
{
	unsigned int median_idx;
	unsigned int n_elements = node.leaf_size();

	//unsigned int axis = node.axis();

	if (n_elements > max_median_sample_size) {
		// random sample
		coors_.resize(max_median_sample_size);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = node.leaf_begin() + ( rand() << 15 | rand() ) % n_elements;
			coors_[i] = elements_[ list_element_index_[ median_idx ] ].projection_center(axis);
		}

	} else {
		// all elements
		coors_.resize(n_elements);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = node.leaf_begin() + i;
			coors_[i] = elements_[ list_element_index_[ median_idx ] ].projection_center(axis);
		}

	}

	unsigned int median_position = (unsigned int)(coors_.size() / 2);
	std::nth_element(coors_.begin(), coors_.begin()+median_position, coors_.end());

	node.median_ = coors_[median_position];
}


void BIHTree::sort_elements(unsigned int &bound1, unsigned int &bound2) {
    unsigned int bound3;
    BIHNode & actual_node = nodes_[queue_.front()];

    bound1 = actual_node.child_[0];
	bound2 = actual_node.child_[0];
	bound3 = actual_node.child_[1];
	while (bound2 != bound3) {
		if (elements_[ list_element_index_[bound2] ].min()(actual_node.axis()) < actual_node.median_) {
			if (elements_[ list_element_index_[bound2] ].max()(actual_node.axis()) > actual_node.median_) {
			    // median in bounding box (element in both ranges)
				bound2++;
			} else {
			    // median after bounding box (element in left range)
				if (bound1 != bound2) {
					std::swap(list_element_index_[bound1], list_element_index_[bound2]);
				}
				bound1++;
				bound2++;
			}
		} else {
		    // median before bounding box (element in right range)
			std::swap(list_element_index_[bound2], list_element_index_[bound3-1]);
			bound3--;
		}
	}
}


void BIHTree::distribute_elements(BIHNode &left_child, BIHNode &right_child) {
    unsigned int bound1, bound2;
    BIHNode & actual_node = nodes_[queue_.front()];
    sort_elements(bound1, bound2);

	// distribute elements into subareas
    left_child.child_[0] = list_element_index_next_.size();
	for (unsigned int i=actual_node.child_[0]; i<bound2; i++) {
		list_element_index_next_.push_back(list_element_index_[i]);
	}
	left_child.child_[1] = list_element_index_next_.size();
	right_child.child_[0] = list_element_index_next_.size();
	for (unsigned int i=bound1; i<actual_node.child_[1]; i++) {
		list_element_index_next_.push_back(list_element_index_[i]);
	}
	right_child.child_[1] = list_element_index_next_.size();
}


void BIHTree::put_leaf_elements() {
	unsigned int lower_bound = in_leaves_.size();
	BIHNode & actual_node = nodes_[queue_.front()];

	for (unsigned int i=actual_node.child_[0]; i<actual_node.child_[1]; i++) {
		in_leaves_.push_back(list_element_index_[i]);
	}

	test_new_level();
	actual_node.child_[0] = lower_bound;
	actual_node.child_[1] = in_leaves_.size();
	queue_.pop_front();
	box_queue_.pop_front();
}


void BIHTree::free_memory() {
	std::deque<unsigned int>().swap(queue_);
	std::deque<BoundingBox>().swap(box_queue_);
	std::vector<unsigned int>().swap(list_element_index_);
	std::vector<unsigned int>().swap(list_element_index_next_);
	std::vector<double>().swap(coors_);
}


void BIHTree::test_new_level() {
	if (nodes_[queue_.front()].child_[1] == list_element_index_.size()) {
		//printf("test splnen\n");
		list_element_index_.swap(list_element_index_next_);
		list_element_index_next_.clear();
	}
}


unsigned int BIHTree::get_element_count() {
	return elements_.size();
}


void BIHTree::find_bounding_box(const BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements)
{
	std::vector<unsigned int>::iterator it;
	searchedElements.clear();
	queue_.clear();

	ASSERT(nodes_.size(), "BIH Tree not created.");

	queue_.push_back(0);
	while (queue_.size()) {
		const BIHNode & node = nodes_[queue_.front()];
		if (node.is_leaf()) {

			DBGMSG( "leaf: %d size %d\n", queue_.front(), node.leaf_size() );
			//START_TIMER("leaf");
			for (unsigned int i=node.leaf_begin(); i<node.leaf_end(); i++) {
				//DBGCOUT( << "in leaf: " << elements_[ in_leaves_[i] ] <<endl );
				if (elements_[ in_leaves_[i] ].intersect(boundingBox)) {
					searchedElements.push_back(in_leaves_[i]);
				}
			}
			//END_TIMER("leaf");
		} else {
			//START_TIMER("recursion");
			if ( boundingBox.min()(node.axis()) < node.median() ) {
				DBGCOUT( << "child 0: " << endl);
				queue_.push_back( node.child(0) );
			}
			if ( boundingBox.max()( node.axis()) > node.median() ) {
				queue_.push_back( node.child(1) );
			}
			//END_TIMER("recursion");
		}
		queue_.pop_front();
	}



	sort(searchedElements.begin(), searchedElements.end());
	it = unique(searchedElements.begin(), searchedElements.end());
	searchedElements.resize( it - searchedElements.begin() );
}


void BIHTree::find_point(const Space<3>::Point &point, std::vector<unsigned int> &searchedElements) {
	unsigned int node_index = 0; // index of actual walking node

	searchedElements.clear();

	while (!nodes_[node_index].is_leaf()) {
		if ( point(nodes_[node_index].axis()) < nodes_[node_index].median() ) {
			node_index = nodes_[node_index].child(0);
		} else {
			node_index = nodes_[node_index].child(1);
		}
	}

	for (unsigned int i=nodes_[node_index].leaf_begin(); i<nodes_[node_index].leaf_end(); i++) {
		if (elements_[ in_leaves_[i] ].contains_point(point)) {
			searchedElements.push_back(in_leaves_[i]);
		}
	}
}


void BIHTree::element_boxes() {
	elements_.resize(mesh_->element.size());
	unsigned int i=0;
	FOR_ELEMENTS(mesh_, element) {
		elements_[i] = element->bounding_box();
	}
}



void BIHTree::get_tree_params(unsigned int &maxDepth, unsigned int &minDepth, double &avgDepth, unsigned int &leafNodesCount,
		unsigned int &innerNodesCount, unsigned int &sumElements) {
	unsigned int sumDepth = 0;
	maxDepth = 0;
	minDepth = 32767;
	leafNodesCount = 0;
	innerNodesCount = 1;

	for (unsigned int i=0; i<nodes_.size(); i++) {
		if (nodes_[i].is_leaf()) {
			if (nodes_[i].depth() > maxDepth) maxDepth = nodes_[i].depth();
			if (nodes_[i].depth() < minDepth) minDepth = nodes_[i].depth();
			sumDepth += nodes_[i].depth();
			++leafNodesCount;
		} else {
			++innerNodesCount;
		}
	}

	avgDepth = (double) sumDepth / (double) leafNodesCount;
	sumElements = in_leaves_.size();
}
