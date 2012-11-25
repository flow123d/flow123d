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

#include "new_mesh/bih_tree.hh"
#include "new_mesh/bih_node.hh"
#include <ctime>

#define DEBUG

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
	unsigned int elementCount, medianCount, medianPosition;
	unsigned int bound1, bound2, bound3, swapIndex, medianIndex;
	unsigned char depth;
	double maxDiff;
	std::vector<double> coors;
	std::vector<unsigned int> listElementIndex;
	std::vector<unsigned int> listElementIndexNext;

	// arma::vec6 stores minimal and maximal coordinations of area
	arma::vec6 areaCoors;

	// temporary vector keeps coordinations of elements stored in queue_
	std::deque<arma::vec6> queueCoors;

	// create root node of tree
	BIHNode bihNode(0);
	bihNode.child_[0] = 0;
	bihNode.child_[1] = mesh_->n_elements();
	nodes_.push_back(bihNode);

	// put indexes of all elements to vector (for first level of tree)
	listElementIndex.resize(mesh_->n_elements());
	for (int i=0; i<mesh_->n_elements(); i++) {
		listElementIndex[i] = i;
	}

	// find minimal and maximal coordination of whole mesh
	Node* node = mesh_->node_vector.begin();
	for (int i=0; i<dimension; i++) {
		areaCoors(i) = node->point()(i);
		areaCoors(i + dimension) = node->point()(i);
	}
	FOR_NODES(mesh_, node ) {
		for (int i=0; i<dimension; i++) {
			areaCoors(i) = std::min( areaCoors(i), node->point()(i) );
			areaCoors(i + dimension) = std::max( areaCoors(i + dimension), node->point()(i) );
		}
	}

	// put index of root node and its coordinations to vectors
	queue_.clear();
	queue_.push_back(0);
	queueCoors.push_back(areaCoors);

	// create tree
	depth = 0;
	while (queue_.size()) {
		BIHNode child[BIHNode::child_count];
		arma::vec6 childCoors[BIHNode::child_count];

		// switch vectors for calculation new level of tree
		if (depth != nodes_[queue_.front()].axes_ - dimension) {
			listElementIndex.swap(listElementIndexNext);
			listElementIndexNext.clear();
			depth = nodes_[queue_.front()].axes_ - dimension;
		}

		// mark node as leaf
		elementCount = nodes_[queue_.front()].get_element_count();
		if (elementCount <= areaElementLimit) {
			put_leaf_elements(listElementIndex);
			queue_.pop_front();
			queueCoors.pop_front();
			continue;
		}

		// Set axes_ class member (select maximal dimension)
		maxDiff = queueCoors.front()(dimension) - queueCoors.front()(0);
		nodes_[queue_.front()].axes_ = 0;
		for (int i=1; i<dimension; i++) {
			if (queueCoors.front()(i + dimension) - queueCoors.front()(i) > maxDiff) {
				maxDiff = queueCoors.front()(i + dimension) - queueCoors.front()(i);
				nodes_[queue_.front()].axes_ = i;
			}
		}

		//select adepts at median
		medianCount = (elementCount >= max_median_count) ? max_median_count : ((elementCount % 2) ? elementCount : elementCount - 1);
		medianPosition = (int)(medianCount/2);
		coors.resize(medianCount);
		for (unsigned int i=0; i<medianCount; i++) {
			medianIndex = nodes_[queue_.front()].child_[0] + ( rand() << 15 | rand() ) % elementCount;
			coors[i] = elements_[ listElementIndex[ medianIndex ] ].get_center()(nodes_[queue_.front()].axes_);
		}

		//select median of the adepts
		std::nth_element(coors.begin(), coors.begin()+medianPosition, coors.end());
		nodes_[queue_.front()].median_ = coors[medianPosition];

		//calculate bounding boxes of subareas and create them
		for (int i=0; i<BIHNode::child_count; i++) {
			for (int j=0; j<dimension; j++) {
				childCoors[i](j) = (j==nodes_[queue_.front()].axes_ && i==1) ? nodes_[queue_.front()].median_
																		: queueCoors.front()(j);
				childCoors[i](j + dimension) = (j==nodes_[queue_.front()].axes_ && i==0) ? nodes_[queue_.front()].median_
																					: queueCoors.front()(j + dimension);
			}

			child[i].set_values(depth+1);
		}

		// sort elements to 3 groups (contained only in left child, contained in both children, contained only in right child)
		bound1 = nodes_[queue_.front()].child_[0];
		bound2 = nodes_[queue_.front()].child_[0];
		bound3 = nodes_[queue_.front()].child_[1];
		unsigned int x1=0, x2=0, x3=0;
		while (bound2 != bound3) {
			if (elements_[ listElementIndex[bound2] ].get_min()(nodes_[queue_.front()].axes_) < nodes_[queue_.front()].median_) {
				if (elements_[ listElementIndex[bound2] ].get_max()(nodes_[queue_.front()].axes_) > nodes_[queue_.front()].median_) {
					bound2++;
					x2++;
				} else {
					if (bound1 != bound2) {
						swapIndex = listElementIndex[bound2];
						listElementIndex[bound2] = listElementIndex[bound1];
						listElementIndex[bound1] = swapIndex;
					}
					bound1++;
					bound2++;
					x1++;
				}
			} else {
				swapIndex = listElementIndex[bound2];
				listElementIndex[bound2] = listElementIndex[bound3-1];
				listElementIndex[bound3-1] = swapIndex;
				bound3--;
				x3++;
			}
		}

		// distribute elements into subareas
		child[0].child_[0] = listElementIndexNext.size();
		for (unsigned int i=nodes_[queue_.front()].child_[0]; i<bound2; i++) {
			listElementIndexNext.push_back(listElementIndex[i]);
		}
		child[0].child_[1] = listElementIndexNext.size();
		child[1].child_[0] = listElementIndexNext.size();
		for (unsigned int i=bound1; i<nodes_[queue_.front()].child_[1]; i++) {
			listElementIndexNext.push_back(listElementIndex[i]);
		}
		child[1].child_[1] = listElementIndexNext.size();

	    //DBGMSG("depth: %d els: %d childs: %d %d %d\n", axes_, get_element_count(),
	    //        n_child_elements, child_[0]->get_element_count(), child_[1]->get_element_count());

		// test count of elements in subareas
		if (child[0].get_element_count() > max_elements_in_child * elementCount ||
		    child[1].get_element_count() > max_elements_in_child * elementCount ||
		    child[0].get_element_count() + child[1].get_element_count() > max_elements_in_children * elementCount) {

			nodes_[queue_.front()].axes_ = depth + dimension;
			put_leaf_elements(listElementIndex);
			queue_.pop_front();
			queueCoors.pop_front();
		    continue;
		}

		// put nodes into vectors
		for (int i=0; i<BIHNode::child_count; i++) {
			nodes_[queue_.front()].child_[i] = nodes_.size();
			queue_.push_back(nodes_.size());
			queueCoors.push_back(childCoors[i]);
			nodes_.push_back(child[i]);
		}

		// remove processed node from queue
		queue_.pop_front();
		queueCoors.pop_front();
	}

	// memory deallocation
	queue_.clear();
	std::deque<unsigned int>().swap(queue_);
}


void BIHTree::put_leaf_elements(std::vector<unsigned int> &listElementId) {
	unsigned int lowerBound = in_leaves_.size();

	for (int i=nodes_[queue_.front()].child_[0]; i<nodes_[queue_.front()].child_[1]; i++) {
		in_leaves_.push_back(listElementId[i]);
	}

	nodes_[queue_.front()].child_[0] = lowerBound;
	nodes_[queue_.front()].child_[1] = in_leaves_.size();
}


unsigned int BIHTree::get_element_count() {
	return elements_.size();
}


void BIHTree::find_elements(BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements)
{
	std::vector<unsigned int>::iterator it;
	searchedElements.clear();
	queue_.clear();

	if (nodes_.size()) {
		queue_.push_back(0);
		while (queue_.size()) {
			if (nodes_[queue_.front()].axes_ >= dimension) {
			    //START_TIMER("leaf");
				for (int i=nodes_[queue_.front()].child_[0]; i<nodes_[queue_.front()].child_[1]; i++) {
					if (elements_[ in_leaves_[i] ].intersection(boundingBox)) {
						searchedElements.push_back(in_leaves_[i]);
					}
				}
				//END_TIMER("leaf");
			} else {
			    //START_TIMER("recursion");
				if ( boundingBox.get_min()(nodes_[ queue_.front() ].axes_) < nodes_[ queue_.front() ].median_ )
					queue_.push_back( nodes_[queue_.front()].child_[0] );
				if ( boundingBox.get_max()(nodes_[ queue_.front() ].axes_) > nodes_[ queue_.front() ].median_ )
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


void BIHTree::element_boxes() {
	elements_.resize(mesh_->element.size());
	FOR_ELEMENTS(mesh_, element) {
		arma::vec3 minCoor = element->node[0]->point();
		arma::vec3 maxCoor = element->node[0]->point();
		int id = element.id();
		for (int i=1; i<element->n_nodes(); i++) {
			Node* node = element->node[i];
			for (int j=0; j<dimension; j++) {
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
		if (nodes_[i].axes_ >= dimension) {
			if (nodes_[i].axes_ - dimension > maxDepth) maxDepth = nodes_[i].axes_ - dimension;
			if (nodes_[i].axes_ - dimension < minDepth) minDepth = nodes_[i].axes_ - dimension;
			sumDepth += nodes_[i].axes_ - dimension;
			++leafNodesCount;
		} else {
			++innerNodesCount;
		}
	}

	avgDepth = (double) sumDepth / (double) leafNodesCount;
	sumElements = in_leaves_.size();
}
