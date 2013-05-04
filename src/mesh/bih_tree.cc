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

const double BIHTree::max_elements_in_child = 0.8;
const double BIHTree::max_elements_in_children = 1.5;


BIHTree::BIHTree(Mesh* mesh, unsigned int areaElementLimit) {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	srand((unsigned)time(0));

	mesh_ = mesh;
	if (areaElementLimit == 0) areaElementLimit = 20;
	nodes_.reserve(2 * mesh_->n_elements() / areaElementLimit);

	//START_TIMER("BIH Tree");

	element_boxes();
	create_tree(areaElementLimit);

	//END_TIMER("BIH Tree");

	xprintf(Msg, " - Tree created\n");
}



BIHTree::~BIHTree() {

}


void BIHTree::create_tree(unsigned int areaElementLimit) {
	create_root_node();

	// create tree
	while (queue_.size()) {
		BIHNode child[BIHNode::child_count];
		arma::vec6 childCoors[BIHNode::child_count];
		ASSERT(queue_.front() < nodes_.size(), "Idx %d out of nodes_ of size %d.\n", queue_.front(), nodes_.size());
		BIHNode & actual_node = nodes_[queue_.front()];

		// mark node as leaf
		if (actual_node.get_element_count() <= areaElementLimit) {
			put_leaf_elements();
			continue;
		}

		unsigned char depth = actual_node.depth();
		set_axes();
		set_median();

		//calculate bounding boxes of subareas and create them
		for (unsigned int i=0; i<BIHNode::child_count; i++) {
			for (unsigned int j=0; j<dimension; j++) {
				childCoors[i](j) = (j==actual_node.axes() && i==1) ? actual_node.median_
																   : queue_coors_.front()(j);
				childCoors[i](j + dimension) = (j==actual_node.axes() && i==0) ? actual_node.median_
																			   : queue_coors_.front()(j + dimension);
			}

			child[i].set_depth(depth+1);
		}

	    distribute_elements(child[0], child[1]);

		// test count of elements in subareas
		if (child[0].get_element_count() > max_elements_in_child * actual_node.get_element_count() ||
		    child[1].get_element_count() > max_elements_in_child * actual_node.get_element_count() ||
		    child[0].get_element_count() + child[1].get_element_count() > max_elements_in_children * actual_node.get_element_count()) {

			list_element_index_next_.erase(list_element_index_next_.begin() + child[0].child_[0],
										   list_element_index_next_.begin() + child[1].child_[1]);
			actual_node.set_depth(depth);
			put_leaf_elements();
		    continue;
		}

		// put nodes into vectors
		test_new_level();
		for (unsigned int i=0; i<BIHNode::child_count; i++) {
		    nodes_[queue_.front()].child_[i] = nodes_.size(); // can not use actual_node, since it can be invalid due to reallocation of nodes_
			queue_.push_back(nodes_.size());
			queue_coors_.push_back(childCoors[i]);
			nodes_.push_back(child[i]);
		}

		// remove processed node from queue
		queue_.pop_front();
		queue_coors_.pop_front();
	}

	free_memory();
}


void BIHTree::create_root_node() {
	// arma::vec6 stores minimal and maximal coordinations of area
	// Mimics BoundingBox functionality.
	// TODO: Possibly add suitable methods to BoundingBox in order to use it here.
	arma::vec6 area_coors;

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
	for (unsigned int i=0; i<dimension; i++) {
		area_coors(i) = node->point()(i);
		area_coors(i + dimension) = node->point()(i);
	}
	FOR_NODES(mesh_, node ) {
		for (unsigned int i=0; i<dimension; i++) {
			area_coors(i) = std::min( area_coors(i), node->point()(i) );
			area_coors(i + dimension) = std::max( area_coors(i + dimension), node->point()(i) );
		}
	}

	// put index of root node and its coordinations to vectors
	queue_.clear();
	queue_.push_back(0);
	queue_coors_.push_back(area_coors);
}


void BIHTree::set_axes() {
	double maxDiff = queue_coors_.front()(dimension) - queue_coors_.front()(0);
	nodes_[queue_.front()].axes_ = 0;
	for (unsigned int i=1; i<dimension; i++) {
		if (queue_coors_.front()(i + dimension) - queue_coors_.front()(i) > maxDiff) {
			maxDiff = queue_coors_.front()(i + dimension) - queue_coors_.front()(i);
			nodes_[queue_.front()].axes_ = i;
		}
	}
}


void BIHTree::set_median() {
	BIHNode & actual_node = nodes_[queue_.front()];

	//select adepts at median
	unsigned int medianIndex;
	unsigned int elementCount = actual_node.get_element_count();
	unsigned int medianCount = (elementCount >= max_median_count) ? max_median_count : ((elementCount % 2) ? elementCount : elementCount - 1);
	unsigned int medianPosition = (unsigned int)(medianCount/2);
	coors_.resize(medianCount);
	for (unsigned int i=0; i<medianCount; i++) {
		medianIndex = actual_node.child_[0] + ( rand() << 15 | rand() ) % elementCount;
		coors_[i] = elements_[ list_element_index_[ medianIndex ] ].get_center()(actual_node.axes());
	}

	//select median of the adepts
	std::nth_element(coors_.begin(), coors_.begin()+medianPosition, coors_.end());
	actual_node.median_ = coors_[medianPosition];
}


void BIHTree::sort_elements(unsigned int &bound1, unsigned int &bound2) {
    unsigned int bound3;
    BIHNode & actual_node = nodes_[queue_.front()];

    bound1 = actual_node.child_[0];
	bound2 = actual_node.child_[0];
	bound3 = actual_node.child_[1];
	while (bound2 != bound3) {
		if (elements_[ list_element_index_[bound2] ].get_min()(actual_node.axes()) < actual_node.median_) {
			if (elements_[ list_element_index_[bound2] ].get_max()(actual_node.axes()) > actual_node.median_) {
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
	queue_coors_.pop_front();
}


void BIHTree::free_memory() {
	std::deque<unsigned int>().swap(queue_);
	std::deque<arma::vec6>().swap(queue_coors_);
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


void BIHTree::find_bounding_box(BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements)
{
	std::vector<unsigned int>::iterator it;
	searchedElements.clear();
	queue_.clear();

	if (nodes_.size()) {
		queue_.push_back(0);
		while (queue_.size()) {
			if (nodes_[queue_.front()].is_leaf()) {
			    //START_TIMER("leaf");
				for (unsigned int i=nodes_[queue_.front()].child_[0]; i<nodes_[queue_.front()].child_[1]; i++) {
					if (elements_[ in_leaves_[i] ].intersection(boundingBox)) {
						searchedElements.push_back(in_leaves_[i]);
					}
				}
				//END_TIMER("leaf");
			} else {
			    //START_TIMER("recursion");
				if ( boundingBox.get_min()(nodes_[ queue_.front() ].axes()) < nodes_[ queue_.front() ].median_ )
					queue_.push_back( nodes_[queue_.front()].child_[0] );
				if ( boundingBox.get_max()(nodes_[ queue_.front() ].axes()) > nodes_[ queue_.front() ].median_ )
					queue_.push_back( nodes_[queue_.front()].child_[1] );
				//END_TIMER("recursion");
			}
			queue_.pop_front();
		}

	}

	sort(searchedElements.begin(), searchedElements.end());
	it = unique(searchedElements.begin(), searchedElements.end());
	searchedElements.resize( it - searchedElements.begin() );
}


void BIHTree::find_point(Point<3> &point, std::vector<unsigned int> &searchedElements) {
	unsigned int node_index = 0; // index of actual walking node

	searchedElements.clear();

	while (!nodes_[node_index].is_leaf()) {
		if ( point(nodes_[node_index].axes()) < nodes_[node_index].median_ ) {
			node_index = nodes_[node_index].child_[0];
		} else {
			node_index = nodes_[node_index].child_[1];
		}
	}

	for (unsigned int i=nodes_[node_index].child_[0]; i<nodes_[node_index].child_[1]; i++) {
		if (elements_[ in_leaves_[i] ].contains_point(point)) {
			searchedElements.push_back(in_leaves_[i]);
		}
	}
}


void BIHTree::element_boxes() {
	elements_.resize(mesh_->element.size());
	FOR_ELEMENTS(mesh_, element) {
		arma::vec3 minCoor = element->node[0]->point();
		arma::vec3 maxCoor = element->node[0]->point();
		for (unsigned int i=1; i<element->n_nodes(); i++) {
			Node* node = element->node[i];
			for (unsigned int j=0; j<dimension; j++) {
				minCoor(j) = std::min(minCoor(j), node->point()(j));
				maxCoor(j) = std::max(maxCoor(j), node->point()(j));
			}
		}
		elements_[element.index()].set_bounds(minCoor, maxCoor);
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
