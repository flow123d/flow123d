/*
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
 * $Id: bih_node.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "new_mesh/bih_node.hh"
#include "new_mesh/bounding_box.hh"
#include "new_mesh/bih_tree.hh"

BIHNode::BIHNode(arma::vec3 minCoordinates, arma::vec3 maxCoordinates, unsigned int depth) {
	//xprintf(Msg, " - BIHNode->BIHNode(arma::vec3, arma::vec3, unsigned int)\n");

	set_values(minCoordinates, maxCoordinates, depth);
}

BIHNode::~BIHNode() {

}


void BIHNode::set_values(arma::vec3 minCoordinates, arma::vec3 maxCoordinates, unsigned int depth) {
	axes_ = 255;
	boundingBox_.set_bounds(minCoordinates, maxCoordinates);
	depth_ = depth;
}


void BIHNode::put_element(unsigned int element_id) {
	element_ids_.push_back(element_id);
}


double BIHNode::get_median_coord(std::vector<BoundingBox> &elements, unsigned int index) {
	unsigned int boundingBoxIndex = element_ids_[index];
	return elements[boundingBoxIndex].get_center()(axes_);
}


void BIHNode::split_distribute(std::vector<BoundingBox> &elements, std::vector<BIHNode> &nodes, unsigned int areaElementLimit) {
	if (get_element_count() <= areaElementLimit) {
		return;
	}

	unsigned int elementCount = get_element_count();
	unsigned int medianCount = (elementCount >= max_median_count) ? max_median_count : ((elementCount % 2) ? elementCount : elementCount - 1);
	unsigned int medianPosition = (int)(medianCount/2);
	double median;
	std:vector<double> coors;
	bool isMaxSplit;
	arma::vec3 diff = boundingBox_.get_max() - boundingBox_.get_min();
	BIHNode child[child_count];

	// Set axes_ class member (select maximal dimension)
	for (int i=0; i<BIHTree::dimension; i++) {
		isMaxSplit = true;
		for (int j=i+1; j<BIHTree::dimension; j++) {
			if (diff(i) < diff(j)) {
				isMaxSplit = false;
				break;
			}
		}
		if (isMaxSplit) {
			axes_ = i;
			break;
		}
	}

	//select adepts at median
	// rand() << 15 + rand();
	coors.resize(medianCount);
	for (unsigned int i=0; i<medianCount; i++) {
		coors[i] = get_median_coord(elements, rand() % elementCount);
	}

	//select median of the adepts
	std::nth_element(coors.begin(), coors.begin()+medianPosition, coors.end());
	median = coors[medianPosition];

	//calculate bounding boxes of subareas and create them
	for (int i=0; i<child_count; i++) {
		arma::vec3 minCoor;
		arma::vec3 maxCoor;
		for (int j=0; j<3; j++) {
			minCoor(j) = (j==axes_ && i==1) ? median : boundingBox_.get_min()(j);
			maxCoor(j) = (j==axes_ && i==0) ? median : boundingBox_.get_max()(j);

		}

		child[i].set_values(minCoor, maxCoor, depth_+1);
	}

	// distribute elements into subareas
    unsigned int n_child_elements=0;
	for (std::vector<unsigned int>::iterator it = element_ids_.begin(); it!=element_ids_.end(); it++) {
		for (int j=0; j<child_count; j++) {
			if (child[j].contains_element(axes_, elements[*it].get_min()(axes_), elements[*it].get_max()(axes_))) {
				child[j].put_element(*it);
				n_child_elements++;
			}

		}
	}
    //DBGMSG("depth: %d els: %d childs: %d %d %d\n", depth_, get_element_count(),
    //        n_child_elements, child_[0]->get_element_count(), child_[1]->get_element_count());

	// test count of elements in subareas
	if (child[0].get_element_count() > 0.8*get_element_count() ||
	    child[1].get_element_count() > 0.8*get_element_count() ||
	    n_child_elements > 1.5 * element_ids_.size()) {

		axes_ = 255;
	    return;
	}

	// put nodes into vector
	child_[0] = nodes.size();
	nodes.push_back(child[0]);
	child_[1] = nodes.size();
	nodes.push_back(child[1]);

	element_ids_.erase(element_ids_.begin(), element_ids_.end());

	nodes[child_[0]].split_distribute(elements, nodes, areaElementLimit);
	nodes[child_[1]].split_distribute(elements, nodes, areaElementLimit);
}


unsigned int BIHNode::get_element_count() {
	return element_ids_.size();
}


bool BIHNode::contains_element(unsigned char coor, double min, double max) {
	return (min < boundingBox_.get_max()(coor)) & (max > boundingBox_.get_min()(coor));
}


void BIHNode::get_tree_params(unsigned int &maxDepth, unsigned int &minDepth, unsigned int &sumDepth, unsigned int &leafNodesCount,
		unsigned int &innerNodesCount, unsigned int &sumElements, std::vector<BIHNode> &nodes) {
	if (axes_ == 255) {
		if (depth_ > maxDepth) maxDepth = depth_;
		if (depth_ < minDepth) minDepth = depth_;
		sumDepth += depth_;
		++leafNodesCount;
		sumElements += element_ids_.size();
	} else {
		++innerNodesCount;
		nodes[child_[0]].get_tree_params(maxDepth, minDepth, sumDepth, leafNodesCount, innerNodesCount, sumElements, nodes);
		nodes[child_[1]].get_tree_params(maxDepth, minDepth, sumDepth, leafNodesCount, innerNodesCount, sumElements, nodes);
	}
}
