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
 * $Id: bounding_interval_hierarchy.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef OCTREE_HH_
#define OCTREE_HH_

#include "system/system.hh"
#include "mesh/mesh.h"
#include "new_mesh/bounding_box.hh"
#include "new_mesh/ngh/include/triangle.h"
#include <armadillo>
#include <algorithm>

class BoundingIntevalHierachy {
public:
    /// count of dimensions
    static const unsigned int dimension = 3;

	/**
	 * Constructor
	 *
	 * Set class members and call functions which create tree
	 * @param mesh Mesh is used for creation the tree
	 */
    BoundingIntevalHierachy(Mesh* mesh);

	/**
	 * Gets suspect elements which can contain point
	 *
	 * @param point Point which is tested if is contained in elements
	 * @param searchedElements vector of suspect elements
	 * @return Negative value if point is outside the area, zero if point is inside
	 */
    int get_element(arma::vec3 &point, std::vector<BoundingBox *> &searchedElements);

	/**
	 * Gets elements which can have intersection with triangle
	 *
	 * @param triangle Triangle which is tested if has intersection
	 * @param searchedElements vector of suspect elements
	 */
    void find_elements(TTriangle &triangle, std::vector<BoundingBox *> &searchedElements);

    /**
     * Add element into elements_ member
     *
     * @param element Added element
     */
	void put_element(BoundingBox* element);

	/**
	 * Get count of elements stored in
	 *
	 * @return Count of elements stored in elements_ member
	 */
    int get_element_count();

private:

    /// limit of elements in area, if count of elements is lesser than value splitting is stopped
    static const unsigned int area_element_limit = 20;
    /// count of subareas - don't change
    static const unsigned int child_count = 2;
    /// count of elements of which is selected median - value must be even
    static const unsigned int area_median_count = 9;

	/**
	 * Constructor
	 *
	 * Set class members and call functions which create children of node
	 * @param mesh Mesh is used for creation the tree
	 * @param minCoordinates Minimal coordinates of BoxElement
	 * @param maxCoordinates Maximal coordinates of BoxElement
	 * @param elementSize Count of elements in parent node
	 * @param splitCoor Coordination of splitting parent area
	 * @param depth Depth of node in tree.
	 */
	BoundingIntevalHierachy(Mesh* mesh, arma::vec3 minCoordinates, arma::vec3 maxCoordinates, int splitCoor, int depth);

    /**
     * Method checks count of elements in area.
     * If count is greater than area_element_limit splits area and distributes elements to subareas.
     */
    void split_distribute();

    /// split area into two subareas by median
    void split_area();
    /// distribute elements into subareas
    void distribute_elements();
    /// create bounding box of area
    void bounding_box();
    /// create bounding boxes of element
    void element_boxes();

    /**
     * Tests if element is contained in area bounding box.
     *
     * @param coor Testing coordination (split coordination of parent area)
     * @param min Minimum value of testing coordination
     * @param max Maximum value of testing coordination
     * @return True if element is contained in area
     */
    bool contains_element(int coor, double min, double max);

    /**
     * Tests if point is contained in area bounding box.
     *
     * @param point Testing point
     * @return True if point is contained in area
     */
    bool contains_point(arma::vec3 &point);

    /// mesh
    Mesh* mesh_;
    /// child nodes
    BoundingIntevalHierachy* child_[child_count];
    /// bounding box of area
    BoundingBox* boundingBox_;
	/// vector of bounding boxes contained in node
    std::vector<BoundingBox *> elements_;
    /// coordination of splitting area
    int splitCoor_;
    /// depth of node
    int depth_;
    /// indicates if node is leaf
    bool leaf_;

};


#endif /* OCTREE_HH_ */
