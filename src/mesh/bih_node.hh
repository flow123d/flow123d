/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    bih_node.hh
 * @brief   
 */

#ifndef BIH_NODE_HH_
#define BIH_NODE_HH_

#include "system/system.hh"
#include "mesh/bounding_box.hh"
#include <armadillo>
#include <algorithm>

class BIHNode {
public:

    /// count of subareas - don't change
    static const unsigned int child_count = 2;
    /// count of dimensions
    static const unsigned char dimension = 3;

    /**
     * Set leaf node.
     */
    void set_leaf(unsigned int begin, unsigned int end, double bound,unsigned int depth) {
    	child_[0]=begin;
    	child_[1]=end;
    	bound_ = bound;
    	axis_=dimension+depth;
    }

    /**
     * Set non-leaf node.
     */
    void set_non_leaf(unsigned int left, unsigned int right, unsigned int axis) {
    	ASSERT(axis < dimension," ");
    	ASSERT(is_leaf(), " "); // first must be leaf node (must set bound_)
    	axis_=axis;
    	child_[0]=left;
    	child_[1]=right;
    }

    /// return true if node is leaf
    bool is_leaf() const
    { return axis_ >= dimension; }

    /// return depth of leaf node

    inline unsigned char depth() const
    {
    	ASSERT( is_leaf(), "Not leaf node.");
    	return axis_ - dimension;
    }

    unsigned int leaf_begin() const
    {
    	ASSERT( is_leaf(), "Not leaf node.");
    	return child_[0];
    }

    unsigned int leaf_end() const
    {
    	ASSERT( is_leaf(), "Not leaf node.");
    	return child_[1];
    }

	/**
	 * Get count of elements stored in
	 *
	 * @return Count of elements contained in node
	 */
    unsigned int leaf_size() const
    {
    	ASSERT( is_leaf(), "Not leaf node.");
    	return child_[1] - child_[0];
    }

    /// return axes (coordination of splitting) of inner node
    unsigned int axis() const
    {
    	ASSERT(!is_leaf(), "Not in branch node.\n");
    	return axis_;
    }

    double bound() const
    {
    	return bound_;
    }

    /// Return index of child node.
    unsigned int child(unsigned int i_child)  const
    {
    	ASSERT(!is_leaf(), "Not in branch node.\n");
    	ASSERT_LESS( i_child, child_count );
    	return child_[i_child];

    }
private:


	/**
	 * Set depth of node to axes_ class members
	 *
	 * @param depth Depth of node in tree.
	 */
    void set_depth(unsigned int depth) {
    	axis_ = depth + dimension;
    }

    /// child nodes indexes
    unsigned int child_[child_count];

    /**
     * A non-leaf node has two childs (left and right). Their bounding boxes are created from the bounding box of parent
     * so that in one direction the left child set max to its bound_ and the right child set its min to bound_.
     * This way we can always repcreate bounding box of every node when traversing the tree from the root.
     */
    double bound_;

    /**
     * Value stores coordination of splitting area for inner nodes or depth for leaf nodes
     *  - values 0,1,2 indicate inner node of tree and coordination of splitting area
     *  - values 3 and greater indicate leaf node of tree and store depth of node (depth = axes_ - 3)
     */
    unsigned char axis_;

};

#endif /* BIH_NODE_HH_ */
