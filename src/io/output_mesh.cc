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
{
    DBGMSG("create outputmesh\n");
    nodes_ = std::make_shared<MeshData<double>>("nodes", OutputDataBase::N_VECTOR);
    connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    offsets_ = std::make_shared<MeshData<unsigned int>>("offsets");
    
    fill_vectors(mesh);
    DBGMSG("create outputmesh\n");
}

OutputMesh::~OutputMesh()
{
}


void OutputMesh::fill_vectors(Mesh* mesh)
{
    const unsigned int n_elements = mesh->n_elements(),
                       n_nodes = mesh->n_nodes();

    nodes_->data_.resize(3*n_nodes);    // suppose 3D coordinates
    nodes_->n_values = n_nodes;
    unsigned int coord_id = 0,  // coordinate id in vector
                 node_id = 0;   // node id
    FOR_NODES(mesh, node) {
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
                 offset = 0,    //offset of node indices of element in node vector
                 li;            //local node index
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            connectivity_->data_.push_back(node->aux);
            connect_id++;
        }
        
        offset += ele->dim() + 1;
        offsets_->data_[ele_id] = offset;
        ele_id++;
        // increase offset by number of nodes of the simplicial element
    }
    connectivity_->data_.shrink_to_fit();
    connectivity_->n_values = connect_id;
}



// void OutputVTK::write_vtk_geometry(void)
// {
//     Mesh *mesh = this->_mesh;
//     ofstream &file = this->_data_file;
// 
//     int tmp;
// 
//     /* Write Points begin*/
//     file << "<Points>" << endl;
//     /* Write DataArray begin */
//     file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
//     /* Write own coordinates */
//     tmp = 0;
//     /* Set floating point precision */
//     file.precision(std::numeric_limits<double>::digits10);
//     FOR_NODES(mesh, node) {
//         node->aux = tmp;   /* store index in the auxiliary variable */
// 
//         file << scientific << node->getX() << " ";
//         file << scientific << node->getY() << " ";
//         file << scientific << node->getZ() << " ";
// 
//         tmp++;
//     }
//     /* Write DataArray end */
//     file << endl << "</DataArray>" << endl;
//     /* Write Points end*/
//     file << "</Points>" << endl;
// }
// 
// void OutputVTK::write_vtk_topology(void)
// {
//     Mesh *mesh = this->_mesh;
//     ofstream &file = this->_data_file;
// 
//     Node* node;
//     //ElementIter ele;
//     unsigned int li;
//     int tmp;
// 
//     /* Write Cells begin*/
//     file << "<Cells>" << endl;
//     /* Write DataArray begin */
//     file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
//     /* Write own coordinates */
//     FOR_ELEMENTS(mesh, ele) {
//         FOR_ELEMENT_NODES(ele, li) {
//             node = ele->node[li];
//             file << node->aux << " ";   /* Write connectivity */
//         }
//     }
//     /* Write DataArray end */
//     file << endl << "</DataArray>" << endl;
// 
//     /* Write DataArray begin */
//     file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
//     /* Write number of nodes for each element */
//     tmp = 0;
//     FOR_ELEMENTS(mesh, ele) {
//         switch(ele->dim()) {
//         case 1:
//             tmp += VTK_LINE_SIZE;
//             break;
//         case 2:
//             tmp += VTK_TRIANGLE_SIZE;
//             break;
//         case 3:
//             tmp += VTK_TETRA_SIZE;
//             break;
//         }
//         file << tmp << " ";
//     }
//     /* Write DataArray end */
//     file << endl << "</DataArray>" << endl;
// 
//     /* Write DataArray begin */
//     file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
//     /* Write type of nodes for each element */
//     FOR_ELEMENTS(mesh, ele) {
//         switch(ele->dim()) {
//         case 1:
//             file << (int)VTK_LINE << " ";
//             break;
//         case 2:
//             file << (int)VTK_TRIANGLE << " ";
//             break;
//         case 3:
//             file << (int)VTK_TETRA << " ";
//             break;
//         }
//     }
//     /* Write DataArray end */
//     file << endl << "</DataArray>" << endl;
// 
//     /* Write Cells end*/
//     file << "</Cells>" << endl;
// }