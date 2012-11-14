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

#include "new_mesh/bih_tree.hh"
#include "new_mesh/bih_node.hh"

#define DEBUG

BIHTree::BIHTree(Mesh* mesh, unsigned int areaElementLimit) {
	xprintf(Msg, " - BIHTree->BIHTree(Mesh, unsigned int)\n");

	mesh_ = mesh;
	if (areaElementLimit == 0) areaElementLimit = 20;
	nodes_.reserve(4 * mesh_->n_elements() / areaElementLimit);

	//START_TIMER("BIH Tree");

	element_boxes();
	root_node(areaElementLimit);

	//END_TIMER("BIH Tree");

	xprintf(Msg, " - Tree created\n");
}



BIHTree::~BIHTree() {

}


void BIHTree::root_node(unsigned int areaElementLimit) {
	Node* node = mesh_->node_vector.begin();
	arma::vec3 point = node->point();

	arma::vec3 minCoordinates = point;
	arma::vec3 maxCoordinates = point;

	FOR_NODES(mesh_, node ) {
		arma::vec3 point = node->point();

		for (int i=0; i<dimension; i++) {
			minCoordinates(i) = std::min(minCoordinates(i), point(i));
			maxCoordinates(i) = std::max(maxCoordinates(i), point(i));
		}
	}

	BIHNode bihNode(minCoordinates, maxCoordinates, 0);
	bihNode.element_ids_.resize(mesh_->n_elements());
	for (int i=0; i<mesh_->n_elements(); i++) {
		bihNode.element_ids_[i] = i;
	}
	nodes_.push_back(bihNode);
	bihNode.split_distribute(elements_, nodes_, areaElementLimit);
}



unsigned int BIHTree::get_element_count() {
	return elements_.size();
}


void BIHTree::find_elements(BoundingBox &boundingBox, std::vector<unsigned int> &searchedElements)
{
	vector<unsigned int>::iterator it;
	searchedElements.clear();
	if (nodes_.size()) {
		nodes_[0].find_elements(boundingBox, searchedElements, elements_, nodes_);
	}

	sort(searchedElements.begin(), searchedElements.end());
	it = unique(searchedElements.begin(), searchedElements.end());
	searchedElements.resize( it - searchedElements.begin() );
}


void BIHTree::element_boxes() {
	elements_.resize(mesh_->element.size());
	FOR_ELEMENTS(mesh_, element) {
		arma::vec3 minCoor = element->node[0]->point();
		arma::vec3 maxCoor = element->node[0]->point();
		int id = element.id();
		for (int i=1; i<element->n_nodes(); i++) {
			Node* node = element->node[i];
			for (int j=0; j<dimension; j++) {
				minCoor(j) = std::min(minCoor(j), node->point()(j));
				maxCoor(j) = std::max(maxCoor(j), node->point()(j));
			}
		}
		elements_[element.index()].set_bounds(minCoor, maxCoor);
		elements_[element.index()].setId(id);
	}
}



void BIHTree::get_tree_params(unsigned int &maxDepth, unsigned int &minDepth, double &avgDepth, unsigned int &leafNodesCount,
		unsigned int &innerNodesCount, unsigned int &sumElements) {
	unsigned int sumDepth = 0;
	maxDepth = 0;
	minDepth = 32767;
	leafNodesCount = 0;
	innerNodesCount = 1;
	sumElements = 0;

	nodes_[0].get_tree_params(maxDepth, minDepth, sumDepth, leafNodesCount, innerNodesCount, sumElements, nodes_);

	avgDepth = (double) sumDepth / (double) leafNodesCount;
}
