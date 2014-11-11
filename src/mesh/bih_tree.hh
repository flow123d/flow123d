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

#include "mesh/bih_node.hh"
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include <armadillo>



/**
 * @brief Class for O(log N) lookup for intersections with a set of bounding boxes.
 *
 * Notes:
 * Assumes spacedim=3. Implementation was designed for arbitrary number of childs per node, but
 * currently it supports max 2 childs per node (binary tree).
 *
 */
class BIHTree {
public:
    /// count of dimensions
    static const unsigned int dimension = 3;
    /// max count of elements to estimate median - value must be even
    static const unsigned int max_median_sample_size = 5;

    /**
	 * Constructor
	 *
	 * Set class members and call functions which create tree
	 * @param mesh  - Mesh used for creation the tree
	 * @param soft_leaf_size_limit - Maximal number of elements stored in a leaf node of BIH tree.
	 */
	BIHTree(Mesh* mesh, unsigned int soft_leaf_size_limit = 20);

	/**
	 * Destructor
	 */
	~BIHTree();

	/**
	 * Get count of elements stored in tree
	 *
	 * @return Count of bounding boxes stored in elements_ member
	 */
    unsigned int get_element_count();

    /**
     * Main bounding box of the whole tree.
     */
    const BoundingBox &tree_box();

	/**
	 * Gets elements which can have intersection with bounding box
	 *
	 * @param boundingBox Bounding box which is tested if has intersection
	 * @param searchedElements vector of ids of suspect elements
	 */
    void find_bounding_box(const BoundingBox &boundingBox, std::vector<unsigned int> &result_list);

	/**
	 * Gets elements which can have intersection with point
	 *
	 * @param point Point which is tested if has intersection
	 * @param searchedElements vector of ids of suspect elements
	 */
    void find_point(const Space<3>::Point &point, std::vector<unsigned int> &result_list);

    /**
     * Get vector of mesh elements bounding boxes
     *
     * @return elements_ vector
     */
    std::vector<BoundingBox> &get_elements() { return elements_; }

protected:
    /// required reduction in size of box to allow further splitting
    static const double size_reduce_factor;

    /// create bounding boxes of element
    void element_boxes();

    /// split tree node given by node_idx, distribute elements to child nodes
    void split_node(const BoundingBox &node_box, unsigned int node_idx);

    /// create child nodes of node given by node_idx
    void make_node(const BoundingBox &box, unsigned int node_idx);

    /**
     * For given node takes projection of centers of bounding boxes of its elements to axis given by
     * @p node::axis()
     * and estimate median of these values. That is optimal split point.
     * Precise median is computed for sets smaller then @p max_median_sample_size
     * estimate from random sample is used for larger sets.
     */
    double estimate_median(unsigned char axis, const BIHNode &node);

    /// mesh
    Mesh* mesh_;
	/// vector of mesh elements bounding boxes
    std::vector<BoundingBox> elements_;
    /// vector of tree nodes
    std::vector<BIHNode> nodes_;
    /// Main bounding box.
    BoundingBox main_box_;
    /// Maximal number of elements stored in a leaf node of BIH tree.
    unsigned int leaf_size_limit;
    /// Maximal count of BIH tree levels
    unsigned int max_n_levels;

    /// vector stored element indexes in leaf nodes
    std::vector<unsigned int> in_leaves_;
    /// temporary vector stored values of coordinations for calculating median
    std::vector<double> coors_;

    // random generator
    std::mt19937	r_gen;

};

#endif /* BIH_TREE_HH_ */
