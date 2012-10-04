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
 * $Id: bih_tree.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef BIH_TREE_HH_
#define BIH_TREE_HH_

#include "new_mesh/bounding_interval_hierarchy.hh"

class BIHTree : public BoundingIntevalHierachy {
	friend class BoundingIntevalHierachy;
public:
	/**
	 * Constructor
	 *
	 * Set class members and call functions which create tree
	 * @param mesh Mesh is used for creation the tree
	 */
	BIHTree(Mesh* mesh);

	/**
	 * Get count of elements stored in
	 *
	 * @return Count of elements stored in element_ids_ member
	 */
    int get_element_count();

	/**
	 * Gets elements which can have intersection with bounding box
	 *
	 * @param boundingBox Bounding box which is tested if has intersection
	 * @param searchedElements vector of ids of suspect elements
	 */
    void find_elements(BoundingBox &boundingBox, std::vector<int> &searchedElements);

protected:
    /// distribute elements into subareas
    void distribute_elements(std::vector<BoundingBox *> elements);
    /// get value of coordination for calculate a median
    double get_median_coord(std::vector<BoundingBox *> elements, int index);
    /// create bounding box of area
    void bounding_box();
    /// create bounding boxes of element
    void element_boxes();

    /// mesh
    Mesh* mesh_;
	/// vector of bounding boxes contained in node
    std::vector<BoundingBox *> elements_;

private:
};

#endif /* BIH_TREE_HH_ */
