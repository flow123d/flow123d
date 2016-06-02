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
 * @file    output_mesh.cc
 * @brief   Classes for auxiliary output mesh.
 */

#include "output_mesh.hh"
#include "mesh/mesh.h"

OutputMesh::OutputMesh(Mesh* mesh)
: orig_mesh_(mesh), discont_data_computed_(false)
{
    DBGMSG("create outputmesh\n");
    nodes_ = std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR);
    connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    offsets_ = std::make_shared<MeshData<unsigned int>>("offsets");
    
    fill_vectors();
    DBGMSG("create outputmesh\n");
}

OutputMesh::~OutputMesh()
{
}


void OutputMesh::fill_vectors()
{
    const unsigned int n_elements = orig_mesh_->n_elements(),
                       n_nodes = orig_mesh_->n_nodes();

    nodes_->data_.resize(3*n_nodes);    // suppose 3D coordinates
    nodes_->n_values = n_nodes;
    unsigned int coord_id = 0,  // coordinate id in vector
                 node_id = 0;   // node id
    FOR_NODES(orig_mesh_, node) {
        node->aux = node_id;   // store node index in the auxiliary variable

        nodes_->data_[coord_id] = node->getX();  coord_id++;
        nodes_->data_[coord_id] = node->getY();  coord_id++;
        nodes_->data_[coord_id] = node->getZ();  coord_id++;
        node_id++;
    }
    
    
    connectivity_->data_.reserve(4*n_elements);  //reserve - suppose all elements being tetrahedra (4 nodes)
    offsets_->data_.resize(n_elements);
    offsets_->n_values = n_elements;
    Node* node;
    unsigned int ele_id = 0,
                 connect_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 li;            // local node index
    FOR_ELEMENTS(orig_mesh_, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            connectivity_->data_.push_back(node->aux);
            connect_id++;
        }
        
        // increase offset by number of nodes of the simplicial element
        offset += ele->dim() + 1;
        offsets_->data_[ele_id] = offset;
        ele_id++;
    }
    connectivity_->data_.shrink_to_fit();
    connectivity_->n_values = connect_id;
}


void OutputMesh::compute_discontinuous_data()
{
    if(discont_data_computed_) return;
    
    discont_nodes_ = std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR);
    discont_connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    
    const unsigned int n_corners = orig_mesh_->n_corners();

    discont_nodes_->data_.resize(3*n_corners);    // suppose 3D coordinates
    discont_nodes_->n_values = n_corners;
 
    discont_connectivity_->data_.resize(n_corners);  //reserve - suppose all elements being tetrahedra (4 nodes)
    discont_connectivity_->n_values = n_corners;
    Node* node;
    unsigned int coord_id = 0,  // coordinate id in vector
                 corner_id = 0, // corner index (discontinous node)
                 li;            // local node index
    FOR_ELEMENTS(orig_mesh_, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            
            discont_nodes_->data_[coord_id] = node->getX();  coord_id++;
            discont_nodes_->data_[coord_id] = node->getY();  coord_id++;
            discont_nodes_->data_[coord_id] = node->getZ();  coord_id++;
            
            discont_connectivity_->data_[corner_id] = corner_id;
            corner_id++;
        }
    }
    
    discont_data_computed_ = true;
}