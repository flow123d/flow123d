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
	friend class BIHNode;
	friend class BIHTree;
public:
    /// count of dimensions
    static const unsigned int dimension = 3;

    /// destructor
    virtual ~BoundingIntevalHierachy();

	/**
	 * Gets suspect elements which can contain point
	 *
	 * @param point Point which is tested if is contained in elements
	 * @param searchedElements vector of suspect elements
	 * @return Negative value if point is outside the area, zero if point is inside
	 */
    int get_element(arma::vec3 &point, std::vector<BoundingBox *> &searchedElements);

	/**
	 * Get count of elements stored in
	 *
	 * @return Count of elements stored in elements_ member
	 */
    virtual int get_element_count() { return -1; }

    /**
     * Get sum of element counts in all leaf nodes of tree
     * Method for gtests
     *
     * @param sum Sums up counts of elements
     */
    virtual void sum_elements_in_leaves(int &sum) {}

    /**
     * Browse tree and get its minimal, maximal and average depth
     * Average depth is counted as sumDepth / leavesCount
     * Method for gtests
     *
     * @param maxDepth Gets maximal depth of tree
     * @param minDepth Gets minimal depth of tree
     * @param sumDepth Gets sum of all depths of tree
     * @param leavesCount Gets count of all leaves of tree
     * @param writeAllDepth Method writes depth in all leaf nodes if value is true
     */
    virtual void get_tree_depth(int &maxDepth, int &minDepth, int &sumDepth, int &leavesCount, bool writeAllDepth) {}

protected:

    /// count of subareas - don't change
    static const unsigned int child_count = 2;
    /// count of elements of which is selected median - value must be even
    static const unsigned int area_median_count = 9;

    /**
     * Empty constructor
     */
    BoundingIntevalHierachy() {}

    /**
     * Method checks count of elements in area.
     * If count is greater than areaElementLimit splits area and distributes elements to subareas.
     */
    void split_distribute(const std::vector<BoundingBox *> &elements, int areaElementLimit);

    /// split area into two subareas by median
    void split_area(const std::vector<BoundingBox *> &elements, int areaElementLimit);
    /// distribute elements into subareas
    virtual void distribute_elements(const std::vector<BoundingBox *> &elements, int areaElementLimit) {}
    /// get value of coordination for calculate a median
    virtual double get_median_coord(const std::vector<BoundingBox *> &elements, int index) { return 0.0; }

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

    /// limit of elements in area, if count of elements is lesser than value splitting is stopped
    //unsigned int area_element_limit_;
    /// child nodes
    BoundingIntevalHierachy* child_[child_count];
    /// bounding box of area
    BoundingBox boundingBox_;
    /// coordination of splitting area
    int splitCoor_;
    /// depth of node
    int depth_;
    /// indicates if node is leaf
    bool leaf_;

};


#endif /* OCTREE_HH_ */
