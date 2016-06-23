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
#include <mesh/ref_element.hh>
#include "fields/field.hh"

#include "output_element.hh"

OutputMesh::OutputMesh(Mesh* mesh, unsigned int max_refinement_level)
:   orig_mesh_(mesh),
    discont_data_computed_(false),
    max_refinement_level_(max_refinement_level)
{
}

OutputMesh::~OutputMesh()
{
}

OutputElementIterator OutputMesh::begin()
{
    return OutputElementIterator(OutputElement(0, this));
}

OutputElementIterator OutputMesh::end()
{
    return OutputElementIterator(OutputElement(offsets_->n_values, this));
}

unsigned int OutputMesh::n_elements()
{
    return offsets_->n_values;
}

unsigned int OutputMesh::n_nodes()
{
    return nodes_->n_values;
}

unsigned int OutputMesh::n_nodes_disc()
{
    //ASSERT_DBG(output_mesh_->discont_data_computed_);
    if(discont_data_computed_)
        return discont_nodes_->n_values;
    else
        return 0;
}



void OutputMesh::create_identical_mesh()
{
//     DBGMSG("create identical outputmesh\n");
    nodes_ = std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR);
    connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    offsets_ = std::make_shared<MeshData<unsigned int>>("offsets");
    
    fill_vectors();
    // if output mesh computed, invalid the discontinuous data
    discont_data_computed_ = false;
}



void OutputMesh::fill_vectors()
{
    const unsigned int n_elements = orig_mesh_->n_elements(),
                       n_nodes = orig_mesh_->n_nodes();

    nodes_->data_.resize(spacedim*n_nodes);    // suppose 3D coordinates
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
    
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_elements);
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
        (*orig_element_indices_)[ele_id] = ele_id;
        ele_id++;
    }
    connectivity_->data_.shrink_to_fit();
    connectivity_->n_values = connect_id;
}


void OutputMesh::compute_discontinuous_data()
{
    ASSERT_DBG(nodes_->n_values > 0);   //continuous data already computed
    if(discont_data_computed_) return;
    
    discont_nodes_ = std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR);
    discont_connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    
    // connectivity = for every element list the nodes => its length corresponds to discontinuous data
    const unsigned int n_corners = connectivity_->n_values;

    discont_nodes_->data_.resize(spacedim*n_corners);
    discont_nodes_->n_values = n_corners;
 
    discont_connectivity_->data_.resize(n_corners);
    discont_connectivity_->n_values = n_corners;

    unsigned int coord_id = 0,  // coordinate id in vector
                 corner_id = 0, // corner index (discontinous node)
                 li;            // local node index

    for(OutputElementIterator it = begin(); it != end(); ++it)
    {
        unsigned int n = it->n_nodes(), 
                     ele_idx = it->idx(),
                     con_off = (*offsets_)[ele_idx];
                     
        for(li = 0; li < n; li++)
        {
            // offset of the first coordinate of the first node of the element in nodes_ vector
            unsigned int off = spacedim * (*connectivity_)[con_off - n + li];
            auto &d = nodes_->data_;
            
            discont_nodes_->data_[coord_id] = d[off];   ++coord_id;
            discont_nodes_->data_[coord_id] = d[off+1]; ++coord_id;
            discont_nodes_->data_[coord_id] = d[off+2]; ++coord_id;
            
            discont_connectivity_->data_[corner_id] = corner_id;
            corner_id++;
        }
    }
    
    discont_data_computed_ = true;
}

        
void OutputMesh::create_refined_mesh(Field<3, FieldValue<3>::Scalar> *error_control_field)
{
    ASSERT(0).error("Not implemented yet.");
    
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>();
    local_nodes_ = std::make_shared<std::vector<double>>();
    
    //FIXME: suggest required capacity
    unsigned int capacity = 16*orig_mesh_->n_elements();
    orig_element_indices_->reserve(capacity);
    local_nodes_->reserve(3*capacity);
    
    // create node indices for connectivity
    unsigned int node_id = 0;   // node id
    FOR_NODES(orig_mesh_, node) {
        node->aux = node_id;   // store node index in the auxiliary variable
        node_id++;
    }
    
    unsigned int coord_id = 0,      // coordinate id in vector
                 connect_id = 0;    // connectivity id
                 
    FOR_ELEMENTS(orig_mesh_, ele) {
        arma::vec3 centre = ele->centre();
        
        AuxElement aux_ele;
        aux_ele.connectivity.resize(ele->n_nodes());
//         aux_ele.coords.resize(ele->n_nodes()*spacedim);
        
        for(unsigned int i=0; i<ele->n_nodes(); i++)
        {
            Node* node = ele->node[i];
            
//             aux_ele.coords[i] = node->getX();
//             aux_ele.coords[i+1] = node->getY();
//             aux_ele.coords[i+2] = node->getZ();
            
            nodes_->data_[coord_id] = node->getX();  coord_id++;
            nodes_->data_[coord_id] = node->getY();  coord_id++;
            nodes_->data_[coord_id] = node->getZ();  coord_id++;
            aux_ele.connectivity[connect_id] = node->aux;     connect_id++;   
        }
        
        //add centre
        nodes_->data_[coord_id] = centre[0];  coord_id++;
        nodes_->data_[coord_id] = centre[1];  coord_id++;
        nodes_->data_[coord_id] = centre[2];  coord_id++;
        aux_ele.connectivity[connect_id] = node_id;    node_id++;
            
        //refinement refinement
        // if(refinement_criterion())
        {
            
            // new element:
            // centre + combination of dim nodes from all nodes
//             refine_aux_element(aux_ele, centre, , ele->dim());
            
        }
    }
    
    orig_element_indices_->shrink_to_fit();
    local_nodes_->shrink_to_fit();
}

template<int dim>
void OutputMesh::refine_aux_element(const OutputMesh::AuxElement& aux_element,
                                    std::vector< OutputMesh::AuxElement >& refinement,
                                    unsigned int& last_node_idx)
{
    static const unsigned int n_subelements = 1 << dim;  //2^dim
    
    // FIXME Use RefElement<>::interact<> from intersections
    static const std::vector<std::vector<unsigned int>> line_nodes[] = {
        {}, //0D
        
        {{0,1}},
        
        {{0,1},
         {0,2},
         {1,2}},
        
        {{0,1},
         {0,2},
         {1,2},
         {0,3},
         {1,3},
         {2,3}}
    };
    
    static const std::vector<unsigned int> conn[] = {
        {}, //0D
        
        {0, 2,
         2, 1},
         
        {0, 3, 4,
         3, 1, 5,
         4, 5, 2,
         3, 5, 4},
         
        {0, 4, 5, 3,
         4, 1, 6, 8,
         5, 6, 2, 9,
         7, 8, 9, 3,
         4, 7, 8, 9,
         4, 6, 8, 9,
         4, 7, 5, 9,
         4, 5, 6, 9}
    }; 
//     DBGMSG("level = %d, %d\n", aux_element.level, max_refinement_level_);
 
    ASSERT_DBG(dim == aux_element.nodes.size()-1);
    
    if( ! refinement_criterion(aux_element)) {
        refinement.push_back(aux_element);
        return;
    }
    
    std::vector<AuxElement> subelements(n_subelements);
    
    // FIXME Use RefElement<>::n_nodes and RefElement<>::n_lines from intersections
    const unsigned int n_old_nodes = dim+1,
                       n_new_nodes = (unsigned int)((dim * (dim + 1)) / 2); // new points are in the center of lines
    
    // auxiliary vectors
    std::vector<Space<spacedim>::Point> nodes = aux_element.nodes;
    std::vector<unsigned int> node_numbering = aux_element.connectivity;
    nodes.reserve(n_old_nodes+n_new_nodes);
    node_numbering.reserve(n_old_nodes+n_new_nodes);

    // create new points in the element
    for(unsigned int e=0; e < n_new_nodes; e++)
    {
        Space<spacedim>::Point p = nodes[line_nodes[dim][e][0]]+nodes[line_nodes[dim][e][1]];
        nodes.push_back( p / 2.0);
//         nodes.back().print();
        
        last_node_idx++;
        node_numbering.push_back(last_node_idx);
    }
   
    
    for(unsigned int i=0; i < n_subelements; i++)
    {
        AuxElement& sub_ele = subelements[i];
        sub_ele.nodes.resize(n_old_nodes);
        sub_ele.connectivity.resize(n_old_nodes);
        sub_ele.level = aux_element.level+1;
        
        // over nodes
        for(unsigned int j=0; j < n_old_nodes; j++)
        {
            unsigned int conn_id = (n_old_nodes)*i + j;
            sub_ele.nodes[j] = nodes[conn[dim][conn_id]];
            sub_ele.connectivity[j] = node_numbering[conn[dim][conn_id]];
        }
        
        refine_aux_element<dim>(sub_ele, refinement, last_node_idx);
    }
}

template void OutputMesh::refine_aux_element<1>(const OutputMesh::AuxElement&,std::vector< OutputMesh::AuxElement >&, unsigned int&);
template void OutputMesh::refine_aux_element<2>(const OutputMesh::AuxElement&,std::vector< OutputMesh::AuxElement >&, unsigned int&);
template void OutputMesh::refine_aux_element<3>(const OutputMesh::AuxElement&,std::vector< OutputMesh::AuxElement >&, unsigned int&);
    
    
    
    
    
//     static
//     std::vector<double> coords = { 0, 0, 0,     // 0
//                                    1, 0, 0,     // 1
//                                    0, 1, 0,     // 2
//                                    0.5, 0, 0,   // 3
//                                    0, 0.5, 0,   // 4
//                                    0.5, 0.5, 0}; // 5
//                                    
//     static
//     std::vector<unsigned int> conn = {0, 3, 4,
//                                       3, 1, 5,
//                                       4, 5, 2,
//                                       3, 5, 4
//     };
    
//     for(unsigned int i=0; i < n_subelements; i++)
//     {
//         AuxElement& sub_ele = subelements[i];
//         
//         // over nodes
//         for(unsigned int j=0; j < 3; j++)   //dim+1
//         {
//             unsigned int conn_id = 3*i + j; //dim+1
//             if(conn[conn_id] < 3) {
//                 sub_ele.connectivity[j] = aux_element.connectivity[conn_id];
//                 sub_ele.coords[j * spacedim] = aux_element.coords[conn_id * spacedim];
//                 sub_ele.coords[j * spacedim +1] = aux_element.coords[conn_id * spacedim +1];
//                 sub_ele.coords[j * spacedim +2] = aux_element.coords[conn_id * spacedim +2];
//             }
//             else {
//                 last_node_idx++;
//                 sub_ele.connectivity[j] = last_node_idx;
//                 sub_ele.coords[j * spacedim] = aux_element.coords[0] + coords[conn_id * spacedim] * (1);
//                 sub_ele.coords[j * spacedim +1] = aux_element.coords[conn_id * spacedim +1];
//                 sub_ele.coords[j * spacedim +2] = aux_element.coords[conn_id * spacedim +2];
//             }
//             
//         }
//     }
//     
//     return subelements;
// }


// std::vector< OutputMesh::AuxElement > OutputMesh::refine_aux_element(OutputMesh::AuxElement& aux_element, 
//                                                                      const Space< spacedim >::Point &centre,
//                                                                      unsigned int centre_idx,
//                                                                      unsigned int dim)
// {
// static const std::vector<std::vector<std::vector<unsigned int>>> side_permutations = 
//                 {   { {0},{1}}, // line
//                     { {0,1}, {0,2}, {1,2}}, //triangle
//                     { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}   //tetrahedron
//                 };
//     
//     unsigned int n_subelements = dim+1;
//     std::vector<AuxElement> subelements(n_subelements);
//     
//     for(unsigned int i=0; i < n_subelements; i++)
//     {
//         AuxElement& sub_ele = subelements[i];
//         
//         for(unsigned int j=0; j < dim; j++)
//         {
//             sub_ele.connectivity[j] = aux_element.connectivity[side_permutations[dim][i][j]];
//             sub_ele.coords[j] = aux_element.coords[side_permutations[dim][i][j]];
//         }
//         sub_ele.connectivity[dim] = centre_idx;
//         sub_ele.coords[dim] = centre;
//     }
//     return subelements;
// }


bool OutputMesh::refinement_criterion(const OutputMesh::AuxElement& ele)
{
    return (ele.level < max_refinement_level_);
}
