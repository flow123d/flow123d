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

BIHTree::BIHTree(Mesh* mesh, unsigned int areaElementLimit) : BoundingIntevalHierachy() {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	mesh_ = mesh;
	if (areaElementLimit == 0) area_element_limit_ = 20;
	else area_element_limit_ = areaElementLimit;
	leaf_ = false;
	depth_ = 0;

	bounding_box();
	element_boxes();
	split_area(elements_);
	distribute_elements(elements_);

	xprintf(Msg, " - Tree created\n");
}



BIHTree::~BIHTree() {

}



void BIHTree::bounding_box() {
	Node* node = mesh_->node_vector.begin();
	arma::vec3 point = node->point();

	arma::vec3 minCoordinates = point;
	arma::vec3 maxCoordinates = point;

	FOR_NODES(mesh_, node ) {
		arma::vec3 point = node->point();

		for (int i=0; i<dimension; i++) {
			minCoordinates(i) = std::min(minCoordinates(i), point(i));
			maxCoordinates(i) = std::max(maxCoordinates(i), point(i));
		}
	}
	boundingBox_.set_bounds(minCoordinates, maxCoordinates);

}



int BIHTree::get_element_count() {
	return elements_.size();
}


void BIHTree::sum_elements_in_leaves(int &sum) {
	sum = 0;
	((BIHNode *)child_[0])->sum_elements_in_leaves(sum);
	((BIHNode *)child_[1])->sum_elements_in_leaves(sum);
}


void BIHTree::distribute_elements(std::vector<BoundingBox *> elements)
{
	int index=0;
	for (std::vector<BoundingBox *>::iterator it = elements_.begin(); it!=elements_.end(); it++) {
		for (int j=0; j<child_count; j++) {
			if (child_[j]->contains_element(splitCoor_, ((BoundingBox*)*it)->get_min()(splitCoor_), ((BoundingBox*)*it)->get_max()(splitCoor_))) {
				((BIHNode *)child_[j])->put_element(index);
			}
		}
		++index;
	}

	((BIHNode *)child_[0])->split_distribute(elements_);
	((BIHNode *)child_[1])->split_distribute(elements_);
}



void BIHTree::find_elements(BoundingBox &boundingBox, std::vector<int> &searchedElements)
{
	searchedElements.clear();
	if (!leaf_) {
		if (child_[0]->boundingBox_.intersection(boundingBox))
		    ((BIHNode *)child_[0])->find_elements(boundingBox, searchedElements, elements_);
		if (child_[1]->boundingBox_.intersection(boundingBox))
		    ((BIHNode *)child_[1])->find_elements(boundingBox, searchedElements, elements_);
	}

	std:sort(searchedElements.begin(), searchedElements.end());

	std::vector<int>::iterator it = searchedElements.begin();
	while (it != (searchedElements.end() - 1)) {
		int idxIter = *it;
		int idxNext = *(it + 1);
		if (idxIter == idxNext) {
			searchedElements.erase(it);
		} else {
			++it;
		}
	}
}



double BIHTree::get_median_coord(std::vector<BoundingBox *> elements, int index) {
	return elements[index]->get_center()(splitCoor_);
}



void BIHTree::element_boxes() {
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
		BoundingBox* boxElement = new BoundingBox(minCoor, maxCoor);
		boxElement->setId(id);
		elements_.push_back(boxElement);
	}
}



void BIHTree::get_tree_depth(int &maxDepth, bool writeAllDepth) {
	maxDepth = 0;
	if (writeAllDepth) xprintf(Msg, " - Depth in leaf nodes:\n");
	((BIHNode *)child_[0])->get_tree_depth(maxDepth, writeAllDepth);
	((BIHNode *)child_[1])->get_tree_depth(maxDepth, writeAllDepth);
	if (writeAllDepth) xprintf(Msg, "\n");
	xprintf(Msg, " - Maximal depth of tree : %d\n", maxDepth);
}
