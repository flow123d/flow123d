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

#include "new_mesh/bih_node.hh"
#include "mesh/mesh.h"
#include <armadillo>

/**
 * @brief Class for O(log N) lookup for intersections with a set of bounding boxes.
 *
 * TODO:
 * - v BIHNode udelat inline metody depth() (osetrit ASSERTEM is_leaf()), is_leaf(), axes() (ASSERT ! is_leaf() )\
 *   ... pouzit create_tree, find_elements
 *
 * - v BIHTree::create_tree - prohozeni vrstev na konec, pokusit se nahradit test hloubky za test konce vrstvy,
 *   mit ASSERT na shodu techto testu, pripadne zjistit pricinu problemu
 *
 *
 * - ROZCLENENI CREATE_tree
 *
 * - more precise documentation
 * - elementy Meshe umi vracet svuj BoundingBox
 *
 *
 */
class BIHTree {
public:
    /// count of dimensions
    static const unsigned int dimension = 3;

    /**
	 * Constructor
	 *
	 * Set class members and call functions which create tree
	 * @param mesh Mesh is used for creation the tree
	 * @param areaElementLimit limit of elements in area
	 */
	BIHTree(Mesh* mesh, unsigned int areaElementLimit = 0);

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
    void find_elements(BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements);

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
    /// max count of elements of which is selected median - value must be even
    static const unsigned int max_median_count = 1023;
    /// value indicates ratio of the number of element in node and number of elements of its child
    static const double max_elements_in_child;
    /// value indicates ratio of the number of element in node and number of elements of its children
    static const double max_elements_in_children;

    /// create bounding boxes of element
    void element_boxes();
    /// create tree
    void create_tree(unsigned int areaElementLimit);
    /// Creates root node of tree, finds its bounding coordinations
    /// and pushes elements to the vector using for creating tree
    void create_root_node();
    /// Set axes_ class member (select maximal dimension) of actually processed node
    void set_axes();
    /// Set median_ class member of actually processed node
    void set_median();
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
	// temporary vector keeps coordinations of elements stored in queue_
	std::deque<arma::vec6> queue_coors_;
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
