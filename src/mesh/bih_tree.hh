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
 * TODO:
 * - rename method find_element, radeji find_bounding_box
 * - new method find_point parametr Space<3>::Point , #include "mesh/point.hh"
 *
 */
class BIHTree {
public:
    /// count of dimensions
    static const unsigned int dimension = 3;
    /// max count of elements to estimate median - value must be even
    static const unsigned int max_median_sample_size = 1023;

    /**
	 * Constructor
	 *
	 * Set class members and call functions which create tree
	 * @param mesh Mesh is used for creation the tree
	 * @param areaElementLimit Maximal number of elements stored in a leaf node of BIH tree.
	 */
	BIHTree(Mesh* mesh, unsigned int areaElementLimit = 20);

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
	 * Gets elements which can have intersection with bounding box
	 *
	 * @param boundingBox Bounding box which is tested if has intersection
	 * @param searchedElements vector of ids of suspect elements
	 */
    void find_bounding_box(const BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements);

	/**
	 * Gets elements which can have intersection with point
	 *
	 * @param point Point which is tested if has intersection
	 * @param searchedElements vector of ids of suspect elements
	 */
    void find_point(const Space<3>::Point &point, std::vector<unsigned int> &searchedElements);

    /**
     * Browse tree and get its typical parameters
     * Method for gtests
     *
     * @param maxDepth Gets maximal depth of tree
     * @param minDepth Gets minimal depth of tree
     * @param avgDepth Gets average depth of tree
     * @param leafNodesCount Gets count of all leaf nodes of tree
     * @param innerNodesCount Gets count of all inner nodes of tree
     * @param elementLeafCount Gets sum of elements contained in all leaf nodes
     */
    void get_tree_params(unsigned int &maxDepth, unsigned int &minDepth, double &avgDepth, unsigned int &leafNodesCount,
    		unsigned int &innerNodesCount, unsigned int &sumElements);

    /**
     * Get vector of mesh elements bounding boxes
     *
     * @return elements_ vector
     */
    std::vector<BoundingBox> &get_elements() { return elements_; }

private:
    /// value indicates ratio of the number of element in node and number of elements of its child
    static const double min_reduce_factor;
    /// value indicates ratio of the number of element in node and number of elements of its children
    static const double max_grow_factor;

    /// create bounding boxes of element
    void element_boxes();
    /**
     * soft_leaf_size_limit - we try to split node if its size (number of elements) is
     * greater then this limit. However we stop branching if number of redundant elements is to big.
     */
    void create_tree(unsigned int soft_leaf_size_limit);
    /// Creates root node of tree, finds its bounding coordinations
    /// and pushes elements to the vector using for creating tree
    void create_root_node();

    /**
     * For given node takes projection of centers of bounding boxes of its elements to axis given by
     * @p node::axis()
     * and estimate median of these values. That is optimal split point.
     * Precise median is computed for sets smaller then @p max_median_sample_size
     * estimate from random sample is used for larger sets.
     */
    void set_node_median(unsigned char axis, BIHNode &node);

    /// Put indexes of elements to in_leaves_ vector if node is marked as leaf
    void put_leaf_elements();
    /// Deallocate memory reserved by vectors
    void free_memory();
    /// Test if processed node is at the end of level during creating of tree
    void test_new_level();
    /**
     * Sort elements in node to 3 groups (contained only in left child, in both children, only in right child)
     *
     * @param bound1 Bound of sorting between first and second group
     * @param bound2 Bound of sorting between second and third group
     */
    void sort_elements(unsigned int &bound1, unsigned int &bound2);
    /**
     * Distribute elements into subareas
     *
     * @param leftChild Left child of actually processed node
     * @param rightChild Right child of actually processed node
     */
    void distribute_elements(BIHNode &left_child, BIHNode &right_child);

    /// mesh
    Mesh* mesh_;
	/// vector of mesh elements bounding boxes
    std::vector<BoundingBox> elements_;
    /// vector of tree nodes
    std::vector<BIHNode> nodes_;


    /// vector stored elements for level-order walk of tree
    std::deque<unsigned int> queue_;


    std::deque<BoundingBox> box_queue_;

    /// vector stored element indexes in leaf nodes
    std::vector<unsigned int> in_leaves_;
    /// temporary vector stored element indexes in last complete level of the BFS tree
    std::vector<unsigned int> list_element_index_;
    /// temporary vector stored element indexes in actually constructed level of the BFS tree
    std::vector<unsigned int> list_element_index_next_;
    /// temporary vector stored values of coordinations for calculating median
    std::vector<double> coors_;

};

#endif /* BIH_TREE_HH_ */
