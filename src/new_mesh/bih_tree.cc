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

BIHTree::BIHTree(Mesh* mesh, unsigned int areaElementLimit) {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	srand((unsigned)time(0));

	mesh_ = mesh;
	if (areaElementLimit == 0) areaElementLimit = 20;
	nodes_.reserve(2 * mesh_->n_elements() / areaElementLimit);

	//START_TIMER("BIH Tree");

	element_boxes();
	root_node(areaElementLimit);

	//END_TIMER("BIH Tree");

	xprintf(Msg, " - Tree created\n");
}



BIHTree::~BIHTree() {

}


void BIHTree::root_node(unsigned int areaElementLimit) {
	BIHNode bihNode(0);
	bihNode.element_ids_.resize(mesh_->n_elements());
	for (int i=0; i<mesh_->n_elements(); i++) {
		bihNode.element_ids_[i] = i;
	}
	nodes_.push_back(bihNode);
	create_tree(areaElementLimit);
}


void BIHTree::create_tree(unsigned int areaElementLimit) {
	unsigned int elementCount, medianCount, medianPosition;
	bool isMaxSplit;
	std:vector<double> coors;

	// arma::vec6 stores minimal and maximal coordinations of area
	arma::vec6 areaCoors;

	// temporary vector keeps coordinations of elements stored in queue_
	std::vector<arma::vec6> queueCoors;

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

	while (queue_.size()) {
		BIHNode child[BIHNode::child_count];
		arma::vec6 childCoors[BIHNode::child_count];
		elementCount = nodes_[queue_[0]].get_element_count();

		if (elementCount <= areaElementLimit) {
			queue_.erase(queue_.begin());
			queueCoors.erase(queueCoors.begin());
			continue;
		}

		medianCount = (elementCount >= max_median_count) ? max_median_count : ((elementCount % 2) ? elementCount : elementCount - 1);
		medianPosition = (int)(medianCount/2);

		// Set axes_ class member (select maximal dimension)
		for (int i=0; i<dimension; i++) {
			isMaxSplit = true;
			for (int j=i+1; j<dimension; j++) {
				if (queueCoors[0](i + dimension) - queueCoors[0](i) < queueCoors[0](j + dimension) - queueCoors[0](j)) {
					isMaxSplit = false;
					break;
				}
			}
			if (isMaxSplit) {
				nodes_[queue_[0]].axes_ = i;
				break;
			}
		}

		//select adepts at median
		coors.resize(medianCount);
		for (unsigned int i=0; i<medianCount; i++) {
			coors[i] = nodes_[queue_[0]].get_median_coord(elements_, ( rand() << 15 | rand() ) % elementCount);
		}

		//select median of the adepts
		std::nth_element(coors.begin(), coors.begin()+medianPosition, coors.end());
		nodes_[queue_[0]].median_ = coors[medianPosition];

		//calculate bounding boxes of subareas and create them
		for (int i=0; i<BIHNode::child_count; i++) {
			arma::vec3 minCoor;
			arma::vec3 maxCoor;
			for (int j=0; j<dimension; j++) {
				minCoor(j) = (j==nodes_[queue_[0]].axes_ && i==1) ? nodes_[queue_[0]].median_ : queueCoors[0](j);
				maxCoor(j) = (j==nodes_[queue_[0]].axes_ && i==0) ? nodes_[queue_[0]].median_ : queueCoors[0](j + dimension);
				childCoors[i](j) = minCoor(j);
				childCoors[i](j + dimension) = maxCoor(j);
			}

			child[i].set_values(nodes_[queue_[0]].depth_+1);
		}

		// distribute elements into subareas
		for (std::vector<unsigned int>::iterator it = nodes_[queue_[0]].element_ids_.begin(); it!=nodes_[queue_[0]].element_ids_.end(); it++) {
			if (elements_[*it].get_min()(nodes_[queue_[0]].axes_) < nodes_[queue_[0]].median_) {
				child[0].put_element(*it);
			}
			if (elements_[*it].get_max()(nodes_[queue_[0]].axes_) > nodes_[queue_[0]].median_) {
				child[1].put_element(*it);
			}
		}
	    //DBGMSG("depth: %d els: %d childs: %d %d %d\n", depth_, get_element_count(),
	    //        n_child_elements, child_[0]->get_element_count(), child_[1]->get_element_count());

		// test count of elements in subareas
		if (child[0].get_element_count() > 0.8*elementCount ||
		    child[1].get_element_count() > 0.8*elementCount ||
		    child[0].get_element_count() + child[1].get_element_count() > 1.5 * nodes_[queue_[0]].element_ids_.size()) {

			nodes_[queue_[0]].axes_ = 255;
			queue_.erase(queue_.begin());
			queueCoors.erase(queueCoors.begin());
		    continue;
		}

		// put nodes into vectors
		for (int i=0; i<BIHNode::child_count; i++) {
			nodes_[queue_[0]].child_[i] = nodes_.size();
			queue_.push_back(nodes_.size());
			queueCoors.push_back(childCoors[i]);
			nodes_.push_back(child[i]);
		}

		nodes_[queue_[0]].element_ids_.erase(nodes_[queue_[0]].element_ids_.begin(), nodes_[queue_[0]].element_ids_.end());
		queue_.erase(queue_.begin());
		queueCoors.erase(queueCoors.begin());
	}
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
			if (nodes_[queue_[0]].axes_ == 255) {
			    //START_TIMER("leaf");
				for (it = nodes_[queue_[0]].element_ids_.begin(); it!=nodes_[queue_[0]].element_ids_.end(); it++) {
					if (elements_[*it].intersection(boundingBox)) {
						searchedElements.push_back(*it);
					}
				}
				//END_TIMER("leaf");
			} else {
			    //START_TIMER("recursion");
				if ( boundingBox.get_min()(nodes_[ queue_[0] ].axes_) < nodes_[ queue_[0] ].median_ )
					queue_.push_back( nodes_[queue_[0]].child_[0] );
				if ( boundingBox.get_max()(nodes_[ queue_[0] ].axes_) > nodes_[ queue_[0] ].median_ )
					queue_.push_back( nodes_[queue_[0]].child_[1] );
				//END_TIMER("recursion");
			}
			queue_.erase(queue_.begin());
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
	sumElements = 0;

	for (unsigned int i=0; i<nodes_.size(); i++) {
		if (nodes_[i].axes_ == 255) {
			if (nodes_[i].depth_ > maxDepth) maxDepth = nodes_[i].depth_;
			if (nodes_[i].depth_ < minDepth) minDepth = nodes_[i].depth_;
			sumDepth += nodes_[i].depth_;
			++leafNodesCount;
			sumElements += nodes_[i].get_element_count();
		} else {
			++innerNodesCount;
		}
	}

	avgDepth = (double) sumDepth / (double) leafNodesCount;
}
