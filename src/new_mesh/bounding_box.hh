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
 * $Id: bounding_box.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef BOX_ELEMENT_HH_
#define BOX_ELEMENT_HH_

#include "system/system.hh"
#include <armadillo>

/**
 * Contains data of bounding box.
 * Used for areas and elements.
 */
class BoundingBox {
public:

	/**
	 * Constructor.
	 *
	 * Set class members
	 * @param minCoor Set value to minCoordinates_
	 * @param maxCoor Set value to maxCoordinates_
	 */
	BoundingBox(arma::vec3 minCoor, arma::vec3 maxCoor);

	/**
	 * Set id of element
	 * NOT USED this method!
	 *
	 * @param id Element id for set
	 */
	void setId(int id);
	/// get id of element
    int getId();
    /// get minimal coordinates of bounding box
    arma::vec3 get_min();
    /// get maximal coordinates of bounding box
    arma::vec3 get_max();
    /// get center coordinates of bounding box
    arma::vec3 get_center();

    /**
     * Detects if box element contains point
     *
     * @param point Testing point
     * @return True if box element contains point
     */
    bool contains_point(arma::vec3 &point);

    bool intersection(BoundingBox &b2);

private:
    /// minimal coordinates of bounding box
    arma::vec3 minCoordinates_;
    /// maximal coordinates of bounding box
    arma::vec3 maxCoordinates_;
    /// id of element
    int elementId_;
};

#endif /* BOX_ELEMENT_HH_ */
