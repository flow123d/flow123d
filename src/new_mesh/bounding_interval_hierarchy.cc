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
#include "new_mesh/ngh/include/tetrahedron.h"
#include "new_mesh/ngh/include/point.h"
#include "new_mesh/ngh/include/intersection.h"
#include <typeinfo>

BoundingIntevalHierachy::~BoundingIntevalHierachy() {
	/*if (boundingBox_ != NULL) {
		delete boundingBox_;
	}*/
	if (!leaf_) {
		delete child_[0];
		delete child_[1];
	}
}



bool BoundingIntevalHierachy::contains_element(int coor, double min, double max) {
	return (min < boundingBox_.get_max()(coor)) & (max > boundingBox_.get_min()(coor));
}



bool BoundingIntevalHierachy::contains_point(arma::vec3 &point) {
	for (int i=0; i<dimension; i++) {
		if ((point(i) < boundingBox_.get_min()(i)) | (point(i) > boundingBox_.get_max()(i))) return false;
	}

	return true;
}



int BoundingIntevalHierachy::get_element(arma::vec3 &point, std::vector<BoundingBox *> &searchedElements) {
	return 0;
	/*if (leaf_) { // NEPOUZIVA SE
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
	}*/
}



void BoundingIntevalHierachy::split_distribute(std::vector<BoundingBox *> elements, int areaElementLimit) {
	if (get_element_count() > areaElementLimit) {
		split_area(elements, areaElementLimit);
		if (!leaf_) distribute_elements(elements, areaElementLimit);
	} else {
		leaf_ = true;
	}
}



void BoundingIntevalHierachy::split_area(std::vector<BoundingBox *> elements, int areaElementLimit) {
	int elementCount = get_element_count();
	int medianPosition = (int)(area_median_count/2);
	double median;
	std:vector<double> coors;
	bool isMaxSplit;
	arma::vec3 diff = boundingBox_.get_max() - boundingBox_.get_min();

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
	coors.resize(area_median_count);
	for (int i=0; i<area_median_count; i++) {
		coors[i] = get_median_coord(elements, rand() % elementCount);
	}

	//select median of the adepts
	std::nth_element(coors.begin(), coors.begin()+medianPosition, coors.end());
	median = coors[medianPosition];

	//calculate bounding boxes of subareas and create them
	if (median == boundingBox_.get_min()(splitCoor_) || median == boundingBox_.get_max()(splitCoor_)) {
		leaf_ = true;
	} else {
		for (int i=0; i<child_count; i++) {
			arma::vec3 minCoor;
			arma::vec3 maxCoor;
			for (int j=0; j<3; j++) {
				minCoor(j) = (j==splitCoor_ && i==1) ? median : boundingBox_.get_min()(j);
				maxCoor(j) = (j==splitCoor_ && i==0) ? median : boundingBox_.get_max()(j);
			}

			child_[i] = new BIHNode(minCoor, maxCoor, splitCoor_, depth_+1, areaElementLimit);
		}
	}
}
