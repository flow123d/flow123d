/*!
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
 * $Id: bounding_interval_hierarchy.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "system/system.hh"
#include "system/sys_vector.hh"
#include "new_mesh/bounding_interval_hierarchy.hh"
#include "new_mesh/bih_node.hh"
#include <typeinfo>

BoundingIntevalHierachy::BoundingIntevalHierachy(Mesh* mesh) {
	xprintf(Msg, " - BoundingIntevalHierachy->BoundingIntevalHierachy(Mesh* mesh)\n");

	mesh_ = mesh;
	leaf_ = false;
	depth_ = 0;

	bounding_box();
	element_boxes();
	split_area();
	distribute_elements();

	xprintf(Msg, " - Tree created\n");
}

BoundingIntevalHierachy::BoundingIntevalHierachy(Mesh* mesh, arma::vec3 minCoordinates, arma::vec3 maxCoordinates, int splitCoor, int depth) {
	//xprintf(Msg, " - BoundingIntevalHierachy->BoundingIntevalHierachy(Mesh*, arma::vec3, arma::vec3, int, int)\n");

	mesh_ = mesh;
	leaf_ = false;
	boundingBox_ = new BoundingBox(minCoordinates, maxCoordinates);
	splitCoor_ = splitCoor;
	depth_ = depth;
}

void BoundingIntevalHierachy::bounding_box() {
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
	boundingBox_ = new BoundingBox(minCoordinates, maxCoordinates);

}

void BoundingIntevalHierachy::element_boxes() {
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
		BoundingBox* boxElement = new BoundingBox(id, minCoor, maxCoor);
		elements_.push_back(boxElement);
	}
}

bool BoundingIntevalHierachy::contains_element(int coor, double min, double max) {
	return (min < boundingBox_->get_max()(coor)) & (max > boundingBox_->get_min()(coor));
}

bool BoundingIntevalHierachy::contains_point(arma::vec3 &point) {
	for (int i=0; i<dimension; i++) {
		if ((point(i) < boundingBox_->get_min()(i)) | (point(i) > boundingBox_->get_max()(i))) return false;
	}

	return true;
}

int BoundingIntevalHierachy::get_element(arma::vec3 &point, std::vector<BoundingBox *> &searchedElements) {
	if (leaf_) {
		for (std::vector<BoundingBox *>::iterator tmp = elements_.begin(); tmp!=elements_.end(); tmp++)
		{
			BoundingBox* boundingBox = *tmp;

			if (boundingBox->contains_point(point)) searchedElements.push_back(boundingBox);
		}
		return 0;
	} else {
		if (child_[0]->contains_point(point)) return child_[0]->get_element(point, searchedElements);
		else if (child_[1]->contains_point(point)) return child_[1]->get_element(point, searchedElements);
		else return -1;
	}
}

void BoundingIntevalHierachy::put_element(BoundingBox* element) {
	elements_.push_back(element);
}

int BoundingIntevalHierachy::get_element_count() {
	return elements_.size();
}

void BoundingIntevalHierachy::split_distribute() {
	if (elements_.size()>area_element_limit) {
		split_area();
		if (!leaf_) distribute_elements();
	} else {
		leaf_ = true;
	}
}

void BoundingIntevalHierachy::split_area() {
	int medianStep = elements_.size() / area_median_count;
	int medianPosition = (int)(area_median_count/2);
	double median;
	double coors[area_median_count];
	bool isMaxSplit;
	arma::vec3 diff = boundingBox_->get_max() - boundingBox_->get_min();

	// Set splitCoor_ class member (select maximal dimension)
	for (int i=0; i<dimension; i++) {
		isMaxSplit = true;
		for (int j=i+1; j<dimension; j++) {
			if (diff(i) < diff(j)) {
				isMaxSplit = false;
				break;
			}
		}
		if (isMaxSplit) {
			splitCoor_ = i;
			break;
		}
	}

	//select adepts at median
	for (int i=0; i<area_median_count; i++) coors[i] = elements_[i * medianStep]->get_center()(splitCoor_);

	//select median of the adepts
	for (int i=0; i<medianPosition; i++) {
		int minIndex=i, maxIndex=i;
		double min=coors[i], max=coors[i], change;

		for (int j=i+1; j<area_median_count-i; j++) {
			if (coors[j]<min) {
				min=coors[j];
				minIndex=j;
			} else if (coors[j]>max) {
				max=coors[j];
				maxIndex=j;
			}
		}
		change = coors[i];
		coors[i] = coors[minIndex];
		coors[minIndex] = change;
		if (maxIndex==i) maxIndex=minIndex;
		change = coors[area_median_count-i-1];
		coors[area_median_count-i-1] = coors[maxIndex];
		coors[maxIndex] = change;
	}

	median = coors[medianPosition];

	//calculate bounding boxes of subareas and create them
	if (median == boundingBox_->get_min()(splitCoor_) || median == boundingBox_->get_max()(splitCoor_)) {
		leaf_ = true;
	} else {
		for (int i=0; i<child_count; i++) {
			arma::vec3 minCoor;
			arma::vec3 maxCoor;
			for (int j=0; j<3; j++) {
				minCoor(j) = (j==splitCoor_ && i==1) ? median : boundingBox_->get_min()(j);
				maxCoor(j) = (j==splitCoor_ && i==0) ? median : boundingBox_->get_max()(j);
			}

			child_[i] = new BoundingIntevalHierachy(mesh_, minCoor, maxCoor, splitCoor_, depth_+1);
		}
	}
}

void BoundingIntevalHierachy::distribute_elements()
{
	for (std::vector<BoundingBox *>::iterator it = elements_.begin(); it!=elements_.end(); it++) {
		BoundingBox* boundingBox = *it;
		for (int j=0; j<child_count; j++) {
			if (child_[j]->contains_element(splitCoor_, boundingBox->get_min()(splitCoor_), boundingBox->get_max()(splitCoor_))) {
				child_[j]->put_element(boundingBox);
			}
		}
	}

	elements_.erase(elements_.begin(), elements_.end());

	child_[0]->split_distribute();
	child_[1]->split_distribute();
}
