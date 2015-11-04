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
 * @file    bih_tree.cc
 * @brief   
 */

#include "mesh/bih_tree.hh"
#include "mesh/bih_node.hh"
#include "system/global_defs.h"
#include <ctime>
#include <stack>

/**
 * Minimum reduction of box size to allow
 * splitting of a node during tree creation.
 */
const double BIHTree::size_reduce_factor = 0.8;


BIHTree::BIHTree(Mesh* mesh, unsigned int soft_leaf_size_limit)
: mesh_(mesh), leaf_size_limit(soft_leaf_size_limit), r_gen(123)
{
	ASSERT(mesh != nullptr, " ");
	max_n_levels = 2*log(mesh->n_elements())/log(2);

	nodes_.reserve(2*mesh_->n_elements() / leaf_size_limit);
	element_boxes();

	in_leaves_.resize(elements_.size());
	for(unsigned int i=0; i<in_leaves_.size(); i++) in_leaves_[i] = i;

	// make root node
	nodes_.push_back(BIHNode());
	nodes_.back().set_leaf(0, in_leaves_.size(), 0, 0);
	// make root box
	Node* node = mesh_->node_vector.begin();
	main_box_ = BoundingBox(node->point(), node->point());
	FOR_NODES(mesh_, node ) {
		main_box_.expand( node->point() );
	}

	make_node(main_box_, 0);
}


BIHTree::~BIHTree() {

}


void BIHTree::split_node(const BoundingBox &node_box, unsigned int node_idx) {
	BIHNode &node = nodes_[node_idx];
	ASSERT(node.is_leaf(), " ");
	unsigned int axis = node_box.longest_axis();
	double median = estimate_median(axis, node);

	// split elements in node according to the median
	auto left = in_leaves_.begin() + node.leaf_begin(); // first of unresolved elements in @p in_leaves_
	auto right = in_leaves_.begin() + node.leaf_end()-1; // last of unresolved elements in @p in_leaves_

	double left_bound=node_box.min(axis); // max bound of the left group
	double right_bound=node_box.max(axis); // min bound of the right group

	while (left != right) {
		if  ( elements_[ *left ].projection_center(axis) < median) {
			left_bound = std::max( left_bound, elements_[ *left ].max(axis) );
			++left;
		}
		else {
			while ( left != right
					&&  elements_[ *right ].projection_center(axis) >= median ) {
				right_bound = std::min( right_bound, elements_[ *right ].min(axis) );
				--right;
			}
			std::swap( *left, *right);
		}
	}
	// in any case left==right is now the first element of the right group

	if ( elements_[ *left ].projection_center(axis) < median) {
		left_bound = std::max( left_bound, elements_[ *left ].max(axis) );
		++left;
		++right;
	} else {
		right_bound = std::min( right_bound, elements_[ *right ].min(axis) );
	}

	unsigned int left_begin = node.leaf_begin();
	unsigned int left_end = left - in_leaves_.begin();
	unsigned int right_end = node.leaf_end();
	unsigned int depth = node.depth()+1;
    // create new leaf nodes and possibly call split_node on them
	// can not use node reference anymore
	nodes_.push_back(BIHNode());
	nodes_.back().set_leaf(left_begin, left_end, left_bound, depth);
	nodes_.push_back(BIHNode());
	nodes_.back().set_leaf(left_end, right_end, right_bound, depth);

	nodes_[node_idx].set_non_leaf(nodes_.size()-2, nodes_.size()-1, axis);
}


void BIHTree::make_node(const BoundingBox &box, unsigned int node_idx) {
	// we must refer to the node by index to prevent seg. fault due to nodes_ reallocation

	split_node(box,node_idx);

	{
		BIHNode &node = nodes_[node_idx];
		BIHNode &child = nodes_[ node.child(0) ];
		BoundingBox node_box(box);
		node_box.set_max(node.axis(), child.bound() );
		if (	child.leaf_size() > leaf_size_limit
			&&  child.depth() < max_n_levels
			&&  ( node.axis() != node_box.longest_axis()
			      ||  node_box.size(node_box.longest_axis()) < box.size(node.axis())  * size_reduce_factor )
			)
		{
				make_node(node_box, node.child(0) );
		}
	}

	{
		BIHNode &node = nodes_[node_idx];
		BIHNode &child = nodes_[ node.child(1) ];
		BoundingBox node_box(box);
		node_box.set_min(node.axis(), child.bound() );
		if (	child.leaf_size() > leaf_size_limit
			&&  child.depth() < max_n_levels
			&&  ( node.axis() != node_box.longest_axis()
			      ||  node_box.size(node_box.longest_axis()) < box.size(node.axis())  * size_reduce_factor )
			)
		{
				make_node(node_box, node.child(1) );
		}
	}
}


double BIHTree::estimate_median(unsigned char axis, const BIHNode &node)
{
	unsigned int median_idx;
	unsigned int n_elements = node.leaf_size();

	if (n_elements > max_median_sample_size) {
		// random sample
		std::uniform_int_distribution<unsigned int> distribution(node.leaf_begin(), node.leaf_end()-1);
		coors_.resize(max_median_sample_size);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = distribution(this->r_gen);

			coors_[i] = elements_[ in_leaves_[ median_idx ] ].projection_center(axis);
		}

	} else {
		// all elements
		coors_.resize(n_elements);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = node.leaf_begin() + i;
			coors_[i] = elements_[ in_leaves_[ median_idx ] ].projection_center(axis);
		}

	}

	unsigned int median_position = (unsigned int)(coors_.size() / 2);
	std::nth_element(coors_.begin(), coors_.begin()+median_position, coors_.end());

	return coors_[median_position];
}


unsigned int BIHTree::get_element_count() {
	return elements_.size();
}


const BoundingBox &BIHTree::tree_box() {
	return main_box_;
}


void BIHTree::find_bounding_box(const BoundingBox &box, std::vector<unsigned int> &result_list)
{
	std::stack<unsigned int, std::vector<unsigned int> > node_stack;
	ASSERT_EQUAL(result_list.size() , 0);

	node_stack.push(0);
	while (! node_stack.empty()) {
		const BIHNode &node = nodes_[node_stack.top()];
		//DBGMSG("node: %d\n", node_stack.top() );
		node_stack.pop();


		if (node.is_leaf()) {

			//START_TIMER("leaf");
			for (unsigned int i=node.leaf_begin(); i<node.leaf_end(); i++) {
				if (elements_[ in_leaves_[i] ].intersect(box)) {

					result_list.push_back(in_leaves_[i]);
				}
			}
			//END_TIMER("leaf");
		} else {
			//START_TIMER("recursion");
			if ( ! box.projection_gt( node.axis(), nodes_[node.child(0)].bound() ) ) {
				// box intersects left group
				node_stack.push( node.child(0) );
			}
			if ( ! box.projection_lt( node.axis(), nodes_[node.child(1)].bound() ) ) {
				// box intersects right group
				node_stack.push( node.child(1) );
			}
			//END_TIMER("recursion");
		}
	}


#ifdef DEBUG_ASSERT
	// check uniqueness of element indexes
	sort(result_list.begin(), result_list.end());
	it = unique(result_list.begin(), result_list.end());
	ASSERT_EQUAL(searsearchedElements.size() , it - result_list.begin());
#endif
}


void BIHTree::find_point(const Space<3>::Point &point, std::vector<unsigned int> &result_list) {
	find_bounding_box(BoundingBox(point), result_list);
}


void BIHTree::element_boxes() {
	elements_.resize(mesh_->element.size());
	unsigned int i=0;
	FOR_ELEMENTS(mesh_, element) {
		elements_[i] = element->bounding_box();

		i++;
	}
}
