/*
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
 * $Id: bih_node.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "new_mesh/bih_node.hh"
#include "new_mesh/bounding_box.hh"
#include "new_mesh/bih_tree.hh"

BIHNode::BIHNode(unsigned int depth) {
	//xprintf(Msg, " - BIHNode->BIHNode(unsigned int)\n");

	set_values(depth);
}

BIHNode::~BIHNode() {

}


void BIHNode::set_values(unsigned int depth) {
	axes_ = depth + BIHTree::dimension;
}


void BIHNode::put_element(unsigned int element_id) {
	element_ids_.push_back(element_id);
}


double BIHNode::get_median_coord(std::vector<BoundingBox> &elements, unsigned int index) {
	unsigned int boundingBoxIndex = element_ids_[index];
	return elements[boundingBoxIndex].get_center()(axes_);
}


unsigned int BIHNode::get_element_count() {
	return element_ids_.size();
}
