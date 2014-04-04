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
 * $Id: bih_tree.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
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
/**
 * Maximum grow factor of total number of elements in childs compared to number of elements in
 * the parent node during tree creation.
 */
//const double BIHTree::max_grow_factor = 1.5;

/*
BIHTree::BIHTree(Mesh* mesh, unsigned int soft_leaf_size_limit) {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	srand((unsigned)time(0));

	mesh_ = mesh;
	nodes_.reserve(2 * mesh_->n_elements() / soft_leaf_size_limit);

	//START_TIMER("BIH Tree");

	element_boxes();
	create_tree(soft_leaf_size_limit);

	free_memory();

	//END_TIMER("BIH Tree");

	xprintf(Msg, " - Tree created\n");
}
*/

BIHTree::BIHTree(Mesh* mesh, unsigned int soft_leaf_size_limit)
: r_gen(123), mesh_(mesh), leaf_size_limit(soft_leaf_size_limit)
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




/*
void BIHTree::create_tree(unsigned int soft_leaf_size_limit) {
	create_root_node();

	// create tree
	while (queue_.size()) {
		BIHNode child[BIHNode::child_count];
		BoundingBox child_box[BIHNode::child_count];

		ASSERT(queue_.front() < nodes_.size(), "Idx %d out of nodes_ of size %d.\n", queue_.front(), nodes_.size());

		BIHNode & actual_node = nodes_[queue_.front()];

		ASSERT(actual_node.is_leaf(), "Not leaf: %d\n",queue_.front() );
		unsigned char depth = actual_node.depth();
		unsigned int actual_node_leaf_size = actual_node.leaf_size();

		if (actual_node_leaf_size <= soft_leaf_size_limit) {
			// mark node as definitive leaf
			put_leaf_elements();
			continue;
		}

		unsigned char splitting_axis = box_queue_.front().longest_axis();
		set_node_median(splitting_axis,actual_node);
		actual_node.axis_ = splitting_axis; // actual_node is not leaf anymore

		//DBGMSG("node idx: %d median: %g axis: %d\n", queue_.front(), actual_node.median_, actual_node.axis_);


		//calculate bounding boxes of subareas and create them
		box_queue_.front().split(actual_node.axis(),actual_node.median(),child_box[0], child_box[1]);
		child[0].set_depth(depth+1);
		child[1].set_depth(depth+1);

	    distribute_elements(child[0], child[1]);

		// test count of elements in subareas
		if (child[0].leaf_size() > min_reduce_factor * actual_node_leaf_size ||
		    child[1].leaf_size() > min_reduce_factor * actual_node_leaf_size ||
		    child[0].leaf_size() + child[1].leaf_size() > max_grow_factor * actual_node_leaf_size) {

			list_element_index_next_.erase(list_element_index_next_.begin() + child[0].child_[0],
										   list_element_index_next_.begin() + child[1].child_[1]);
			actual_node.set_depth(depth);
			put_leaf_elements();
		    continue;
		}

		// put nodes into vectors
		test_new_level();
		for (unsigned int i=0; i<BIHNode::child_count; i++) {
			// can not use actual_node, since it can be invalid due to reallocation of nodes_

		    nodes_[queue_.front()].child_[i] = nodes_.size();
			queue_.push_back(nodes_.size());
			box_queue_.push_back(child_box[i]);
			nodes_.push_back(child[i]);
		}

		// remove processed node from queue
		queue_.pop_front();
		box_queue_.pop_front();
	}

}*/


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
			//DBGMSG("left_bound: %g\n", left_bound);
			++left;
		}
		else {
			while ( left != right
					&&  elements_[ *right ].projection_center(axis) >= median ) {
				right_bound = std::min( right_bound, elements_[ *right ].min(axis) );
				//DBGMSG("right_bound: %g\n", right_bound);
				--right;
			}
			std::swap( *left, *right);
		}
	}
	// in any case left==right is now the first element of the right group

	if ( elements_[ *left ].projection_center(axis) < median) {
		left_bound = std::max( left_bound, elements_[ *left ].max(axis) );
		//DBGMSG("left_bound: %g\n", left_bound);
		++left;
		++right;
	} else {
		right_bound = std::min( right_bound, elements_[ *right ].min(axis) );
		//DBGMSG("right_bound: %g\n", right_bound);
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

/*
void BIHTree::call_recursion(BoundingBox parent_box, BIHNode &parent, BoundingBox child_box, BIHNode &child) {

}*/

void BIHTree::make_node(const BoundingBox &box, unsigned int node_idx) {
	// we must refer to the node by index to prevent seg. fault due to nodes_ reallocation

	split_node(box,node_idx);

#ifdef DEBUG_MESSAGES
//	cout << "  branch node, axis: " << node.axis() << " bound: " << node.bound() << endl;
#endif
	{
		BIHNode &node = nodes_[node_idx];
		BIHNode &child = nodes_[ node.child(0) ];
//		DBGMSG("Node: %d\n",node.child(0));
		BoundingBox node_box(box);
		node_box.set_max(node.axis(), child.bound() );
		if (	child.leaf_size() > leaf_size_limit
			&&  child.depth() < max_n_levels
			&&  ( node.axis() != node_box.longest_axis()
			      ||  node_box.size(node_box.longest_axis()) < box.size(node.axis())  * size_reduce_factor )
			)
		{
				make_node(node_box, node.child(0) );
		} else {
#ifdef DEBUG_MESSAGES
//			cout << "  leaf node, depth: " << child.depth() << " leaf range: " << child.leaf_begin() << " " << child.leaf_end() << endl;
//			cout << "  " << node_box << endl;
#endif
		}
	}

	{
		BIHNode &node = nodes_[node_idx];
		BIHNode &child = nodes_[ node.child(1) ];
//		DBGMSG("Node: %d\n",node.child(1));
		BoundingBox node_box(box);
		node_box.set_min(node.axis(), child.bound() );
		if (	child.leaf_size() > leaf_size_limit
			&&  child.depth() < max_n_levels
			&&  ( node.axis() != node_box.longest_axis()
			      ||  node_box.size(node_box.longest_axis()) < box.size(node.axis())  * size_reduce_factor )
			)
		{
				make_node(node_box, node.child(1) );
		} else {
#ifdef DEBUG_MESSAGES
//			cout << "  leaf node, depth: " << child.depth() << " leaf range: " << child.leaf_begin() << " " << child.leaf_end() <<endl;
//			cout << "  " <<node_box;
#endif
		}
	}
}

/*
void BIHTree::create_root_node() {

	BIHNode bih_node(0);
	bih_node.child_[0] = 0;
	bih_node.child_[1] = mesh_->n_elements();
	nodes_.push_back(bih_node);

	// put indexes of all elements to vector (for first level of tree)
	list_element_index_.resize(mesh_->n_elements());
	for (unsigned int i=0; i<mesh_->n_elements(); i++) {
		list_element_index_[i] = i;
	}

	// find minimal and maximal coordination of whole mesh
	Node* node = mesh_->node_vector.begin();
	main_box_ = BoundingBox(node->point(), node->point());

	FOR_NODES(mesh_, node ) {
		main_box_.expand( node->point() );
	}

	// put index of root node and its coordinations to vectors
	queue_.clear();
	queue_.push_back(0);
	box_queue_.push_back(main_box_);
}

*/

/*
void BIHTree::set_node_median(unsigned char axis, BIHNode &node)
{
	unsigned int median_idx;
	unsigned int n_elements = node.leaf_size();

	if (n_elements > max_median_sample_size) {
		// random sample
		coors_.resize(max_median_sample_size);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = node.leaf_begin() + ( rand() << 15 | rand() ) % n_elements;
			coors_[i] = elements_[ list_element_index_[ median_idx ] ].projection_center(axis);
		}

	} else {
		// all elements
		coors_.resize(n_elements);
		for (unsigned int i=0; i<coors_.size(); i++) {
			median_idx = node.leaf_begin() + i;
			coors_[i] = elements_[ list_element_index_[ median_idx ] ].projection_center(axis);
		}

	}

	unsigned int median_position = (unsigned int)(coors_.size() / 2);
	std::nth_element(coors_.begin(), coors_.begin()+median_position, coors_.end());

	node.median_ = coors_[median_position];
}*/

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
			//DBGMSG(" random idx: %d, coor: %g\n", median_idx, coors_[i]);
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

	//DBGMSG("median pos:  %d %g\n", median_position, coors_[median_position]);
	return coors_[median_position];
}
/*
void BIHTree::sort_elements(unsigned int &bound1, unsigned int &bound2) {
    unsigned int bound3;
    BIHNode & actual_node = nodes_[queue_.front()];

    // Four groups:
    // <0, bound1) - left elements
    // <bound1, bound2) - elements in both domains
    // <bound2, bound3) - not yet processed elements
    // <bound3, ..      - right elements
    bound1 = actual_node.child_[0];
	bound2 = actual_node.child_[0];
	bound3 = actual_node.child_[1];
	while (bound2 != bound3) {
		if (elements_[ list_element_index_[bound2] ].min()(actual_node.axis()) < actual_node.median_) {
			if (elements_[ list_element_index_[bound2] ].max()(actual_node.axis()) > actual_node.median_) {
			    // median in bounding box (element in both ranges)
				bound2++;
			} else {
			    // median after bounding box (element in left range)
				if (bound1 != bound2) {
					std::swap(list_element_index_[bound1], list_element_index_[bound2]);
				}
				bound1++;
				bound2++;
			}
		} else {
		    // median before bounding box (element in right range)
			std::swap(list_element_index_[bound2], list_element_index_[bound3-1]);
			bound3--;
		}
	}
}*/

/*
void BIHTree::distribute_elements(BIHNode &left_child, BIHNode &right_child) {
    unsigned int bound1, bound2;
    BIHNode & actual_node = nodes_[queue_.front()];
    sort_elements(bound1, bound2);

	// distribute elements into subareas
    left_child.child_[0] = list_element_index_next_.size();
	for (unsigned int i=actual_node.child_[0]; i<bound2; i++) {
		list_element_index_next_.push_back(list_element_index_[i]);
	}
	left_child.child_[1] = list_element_index_next_.size();
	right_child.child_[0] = list_element_index_next_.size();
	for (unsigned int i=bound1; i<actual_node.child_[1]; i++) {
		list_element_index_next_.push_back(list_element_index_[i]);
	}
	right_child.child_[1] = list_element_index_next_.size();
}


void BIHTree::put_leaf_elements() {
	unsigned int lower_bound = in_leaves_.size();
	BIHNode & actual_node = nodes_[queue_.front()];

	for (unsigned int i=actual_node.child_[0]; i<actual_node.child_[1]; i++) {
		in_leaves_.push_back(list_element_index_[i]);
	}

	test_new_level();
	actual_node.child_[0] = lower_bound;
	actual_node.child_[1] = in_leaves_.size();
	queue_.pop_front();
	box_queue_.pop_front();
}*/

/*
void BIHTree::free_memory() {
	std::deque<unsigned int>().swap(queue_);
	std::deque<BoundingBox>().swap(box_queue_);
	std::vector<unsigned int>().swap(list_element_index_);
	std::vector<unsigned int>().swap(list_element_index_next_);
	std::vector<double>().swap(coors_);
}


void BIHTree::test_new_level() {
	if (nodes_[queue_.front()].child_[1] == list_element_index_.size()) {
		//printf("test splnen\n");
		list_element_index_.swap(list_element_index_next_);
		list_element_index_next_.clear();
	}
}
*/


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

			//DBGMSG( "leaf-size %d\n", node.leaf_size() );
			//START_TIMER("leaf");
			for (unsigned int i=node.leaf_begin(); i<node.leaf_end(); i++) {
				//DBGMSG("check i: %d i_ele: %d id: %d\n", i, in_leaves_[i], mesh_->element( in_leaves_[i] ).id());
				//DBGCOUT( << "in leaf: " << elements_[ in_leaves_[i] ] <<endl );
				if (elements_[ in_leaves_[i] ].intersect(box)) {

					//DBGMSG("el_id: %d\n" , mesh_->element(in_leaves_[i]).id() );
					result_list.push_back(in_leaves_[i]);
				}
			}
			//END_TIMER("leaf");
		} else {
			//DBGMSG("axis: %d bound left: %g right: %g\n",node.axis(),nodes_[node.child(0)].bound(), nodes_[node.child(1)].bound());
			//START_TIMER("recursion");
			if ( ! box.projection_gt( node.axis(), nodes_[node.child(0)].bound() ) ) {
				//DBGMSG("push:%d\n", node.child(0));
				// box intersects left group
				node_stack.push( node.child(0) );
			}
			if ( ! box.projection_lt( node.axis(), nodes_[node.child(1)].bound() ) ) {
				//DBGMSG("push:%d\n", node.child(1));
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

/*

void BIHTree::get_tree_params(unsigned int &maxDepth, unsigned int &minDepth, double &avgDepth, unsigned int &leafNodesCount,
		unsigned int &innerNodesCount, unsigned int &sumElements) {
	unsigned int sumDepth = 0;
	maxDepth = 0;
	minDepth = 32767;
	leafNodesCount = 0;
	innerNodesCount = 1;

	for (unsigned int i=0; i<nodes_.size(); i++) {
		if (nodes_[i].is_leaf()) {
			if (nodes_[i].depth() > maxDepth) maxDepth = nodes_[i].depth();
			if (nodes_[i].depth() < minDepth) minDepth = nodes_[i].depth();
			sumDepth += nodes_[i].depth();
			++leafNodesCount;
		} else {
			++innerNodesCount;
		}
	}

	avgDepth = (double) sumDepth / (double) leafNodesCount;
	sumElements = in_leaves_.size();
}
*/
