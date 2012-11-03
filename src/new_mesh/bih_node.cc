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

BIHNode::BIHNode(arma::vec3 minCoordinates, arma::vec3 maxCoordinates, int splitCoor, int depth, unsigned int areaElementLimit) : BoundingIntevalHierachy() {
	//xprintf(Msg, " - BIHNode->BIHNode(arma::vec3, arma::vec3, int, int, unsigned int)\n");

	leaf_ = true;
	boundingBox_.set_bounds(minCoordinates, maxCoordinates);
	splitCoor_ = splitCoor;
	depth_ = depth;
	child_[0]=NULL;
	child_[1]=NULL;
}

BIHNode::~BIHNode() {

}

void BIHNode::put_element(int element_id) {
	element_ids_.push_back(element_id);
}

double BIHNode::get_median_coord(const std::vector<BoundingBox *> &elements, int index) {
	int boundingBoxIndex = element_ids_[index];
	return elements[boundingBoxIndex]->get_center()(splitCoor_);
}

<<<<<<< .mine
void BIHNode::distribute_elements(const std::vector<BoundingBox *> &elements, int areaElementLimit) {
    unsigned n_child_elements=0;
=======
void BIHNode::distribute_elements(const std::vector<BoundingBox *> &elements, int areaElementLimit) {
>>>>>>> .r1944
	for (std::vector<int>::iterator it = element_ids_.begin(); it!=element_ids_.end(); it++) {
		for (int j=0; j<child_count; j++) {
			if (child_[j]->contains_element(splitCoor_, elements[*it]->get_min()(splitCoor_), elements[*it]->get_max()(splitCoor_))) {
				((BIHNode *)child_[j])->put_element(*it);
				n_child_elements++;
			}

		}
	}
    //DBGMSG("depth: %d els: %d childs: %d %d %d\n", depth_, get_element_count(),
    //        n_child_elements, child_[0]->get_element_count(), child_[1]->get_element_count());

	if (child_[0]->get_element_count() > 0.8*get_element_count() ||
	    child_[1]->get_element_count() > 0.8*get_element_count() ||
	    n_child_elements > 1.5 * element_ids_.size()) {

	    delete child_[0];
	    delete child_[1];
	    return;
	}
	leaf_=false;

	element_ids_.erase(element_ids_.begin(), element_ids_.end());

	((BIHNode *)child_[0])->split_distribute(elements, areaElementLimit);
	((BIHNode *)child_[1])->split_distribute(elements, areaElementLimit);
}

void BIHNode::find_elements(BoundingBox &boundingBox, std::vector<int> &searchedElements,const  std::vector<BoundingBox *> &meshElements) {
	if (leaf_) {
	    //START_TIMER("leaf");
		for (std::vector<int>::iterator it = element_ids_.begin(); it!=element_ids_.end(); it++) {
			if (meshElements[*it]->intersection(boundingBox)) {
				searchedElements.push_back(*it);
			}
		}
		//END_TIMER("leaf");
	} else {
	    //START_TIMER("recursion");
		if (child_[0]->boundingBox_.intersection(boundingBox))
		    ((BIHNode *)child_[0])->find_elements(boundingBox, searchedElements, meshElements);
		if (child_[1]->boundingBox_.intersection(boundingBox))
		    ((BIHNode *)child_[1])->find_elements(boundingBox, searchedElements, meshElements);
		//END_TIMER("recursion");
	}
}

int BIHNode::get_element_count() {
	return element_ids_.size();
}

void BIHNode::sum_elements_in_leaves(int &sum) {
	if (leaf_) {
		sum += element_ids_.size();
	} else {
		((BIHNode *)child_[0])->sum_elements_in_leaves(sum);
		((BIHNode *)child_[1])->sum_elements_in_leaves(sum);
	}
}

void BIHNode::get_tree_depth(int &maxDepth, int &minDepth, int &sumDepth, int &leavesCount, bool writeAllDepth) {
	if (leaf_) {
		if (writeAllDepth) xprintf(Msg, "%d - ", depth_);
		if (depth_ > maxDepth) maxDepth = depth_;
		if (depth_ < minDepth) minDepth = depth_;
		sumDepth += depth_;
		++leavesCount;
	} else {
		((BIHNode *)child_[0])->get_tree_depth(maxDepth, minDepth, sumDepth, leavesCount, writeAllDepth);
		((BIHNode *)child_[1])->get_tree_depth(maxDepth, minDepth, sumDepth, leavesCount, writeAllDepth);
	}
}
