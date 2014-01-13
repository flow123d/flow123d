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
 * $Id: bounding_box.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "system/system.hh"
#include "mesh/bounding_box.hh"


BoundingBox::BoundingBox() {

}

BoundingBox::BoundingBox(arma::vec3 minCoor, arma::vec3 maxCoor) {
	minCoordinates_ = minCoor;
	maxCoordinates_ = maxCoor;
}



BoundingBox::BoundingBox(const vector<arma::vec3> &points) {
	auto it = points.begin();
	maxCoordinates_ = minCoordinates_ = *it;
	for(++it; it != points.end(); ++it) {
		for(unsigned int j=0; j<3; j++) {
			minCoordinates_(j) = std::min( minCoordinates_(j), (*it)[j] );
			maxCoordinates_(j) = std::max( maxCoordinates_(j), (*it)[j] );
		}
	}
}


void BoundingBox::set_bounds(arma::vec3 minCoor, arma::vec3 maxCoor) {
	minCoordinates_ = minCoor;
	maxCoordinates_ = maxCoor;
}

const arma::vec3 BoundingBox::get_min() const {
	return minCoordinates_;
}

const arma::vec3 BoundingBox::get_max() const {
	return maxCoordinates_;
}

arma::vec3 BoundingBox::get_center() const {
	return (maxCoordinates_ + minCoordinates_) / 2;
}

bool BoundingBox::contains_point(const Space<3>::Point &point) const {
	for (unsigned int i=0; i<dimension; i++) {
		if ((point(i) < minCoordinates_(i)) | (point(i) > maxCoordinates_(i))) return false;
	}

	return true;
}

bool BoundingBox::intersection(const BoundingBox &b2) const {
	for (unsigned int i=0; i<dimension; i++) {
		if ((minCoordinates_(i) > b2.get_max()(i)) | (maxCoordinates_(i) < b2.get_min()(i))) return false;
	}
	return true;
}
