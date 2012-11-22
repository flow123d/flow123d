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
 * - v pruchodech do sirky pouzit std::deque
 *   po vytvoreni stromu, provest uvolneni pameti pomoci swap + temporary.
 *   .. podobne pro queueCoors
 *
 * - pro ulozeni hloubky pouzit prostor v axes_
 *
 * - zjednoduseni hledani max. dimenze v create_tree
 *
 * - pro konstanty 0.8 a 1.5 udelat const staticke promenne v BIHTree a vysvetlit jejich vyznam
 *
 *
 * - indexy bounding boxu v listech jsou ulozeny v jednom spolecnem vektoru vector<unsigned int> in_leaves
 *
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
	 * Get count of elements stored in
	 *
	 * @return Count of elements stored in element_ids_ member
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
    /// create root node of tree
    void root_node(unsigned int areaElementLimit);
    /// create tree
    void create_tree(unsigned int areaElementLimit);

    /// mesh
    Mesh* mesh_;
	/// vector of mesh elements bounding boxes
    std::vector<BoundingBox> elements_;
    /// vector of tree nodes
    std::vector<BIHNode> nodes_;
    /// vector stored elements for level-order walk of tree
    std::deque<unsigned int> queue_;

};

#endif /* BIH_TREE_HH_ */
