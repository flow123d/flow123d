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
 * $Id: bih_node.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef BIH_NODE_HH_
#define BIH_NODE_HH_

#include "system/system.hh"
#include "mesh/bounding_box.hh"
#include <armadillo>
#include <algorithm>
#include "input/accessors.hh"

class BIHNode {
	friend class BIHTree;
public:

	/**
	 * Destructor
	 */
	~BIHNode();

	/**
	 * Get count of elements stored in
	 *
	 * @return Count of elements contained in node
	 */
    unsigned int get_element_count();

    /// return true if node is leaf
    inline bool is_leaf()
    	{ return axes_ >= dimension; }

    /// return depth of leaf node
    inline unsigned char depth()
    {
    	ASSERT(is_leaf(), "Call depth() function for inner node.\n");
    	return axes_ - dimension;
    }

    /// return axes (coordination of splitting) of inner node
    inline unsigned char axes()
    {
    	ASSERT(!is_leaf(), "Call axes() function for leaf node.\n");
    	return axes_;
    }

private:
    /// max count of elements of which is selected median - value must be even
    static const unsigned int max_median_count = 1023;
    /// count of subareas - don't change
    static const unsigned int child_count = 2;
    /// count of dimensions
    static const unsigned char dimension = 3;

    /**
     * Empty constructor
     */
    BIHNode() { }

    /**
	 * Constructor
	 *
	 * Set class members
	 * @param depth Depth of node in tree.
	 */
	BIHNode(unsigned int depth);

	/**
	 * Set depth of node to axes_ class members
	 *
	 * @param depth Depth of node in tree.
	 */
	void set_depth(unsigned int depth);

    /// child nodes indexes
    unsigned int child_[child_count];
    /// value of median which splits area to subareas (coordination is getting by axes_)
    double median_;
    /**
     * Value stores coordination of splitting area for inner nodes or depth for leaf nodes
     *  - values 0,1,2 indicate inner node of tree and coordination of splitting area
     *  - values 3 and greater indicate leaf node of tree and store depth of node (depth = axes_ - 3)
     */
    unsigned char axes_;

};

#endif /* BIH_NODE_HH_ */
