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
 * @file    duplicate_nodes.cc
 * @ingroup mesh
 * @brief   Construction of mesh structure with nodes duplicated at fractures.
 */

#include <set>
#include <queue>

#include "system/sys_profiler.hh"
#include "mesh/duplicate_nodes.h"
#include "mesh/mesh.h"
#include "mesh/node_accessor.hh"
#include "mesh/accessors.hh"




MeshObject::MeshObject(unsigned int dim)
    : dim_(dim)
{
  for (unsigned int i=0; i<4; i++)
    faces[i] = nullptr;
}


DuplicateNodes::DuplicateNodes(Mesh* mesh)
: mesh_(mesh),
  n_duplicated_nodes_(0)
{
  init_nodes();
  init_from_edges();
  init_from_elements();
  duplicate_nodes();
}



void DuplicateNodes::init_nodes()
{
  // initialize the nodes
  n_duplicated_nodes_ = mesh_->n_nodes();
  node_dim_.resize(n_duplicated_nodes_);
}



void DuplicateNodes::init_from_edges()
{
  // initialize the objects from edges
  for (auto edge : mesh_->edge_range()) {
    MeshObject obj(edge.side(0)->dim());
    for (unsigned int i=0; i<edge.side(0)->dim()+1; ++i)
      obj.nodes[i] = edge.side(0)->node(i).idx();
    obj_4_edg_.push_back(objects_[edge.side(0)->dim()].size());
    objects_[edge.side(0)->dim()].push_back(obj);
  }
}



void DuplicateNodes::init_from_elements()
{
  // initialize the objects from elements
  for ( auto ele : mesh_->elements_range() ) {
    MeshObject l(ele->dim());
    for (unsigned int i=0; i<ele->dim()+1; i++)
      l.nodes[i] = ele->node_idx(i);
    for (unsigned int i=0; i<ele->dim()+1; i++)
      l.faces[i] = &objects_[ele->dim()-1][obj_4_edg_[ele->edge_idx(i)]];
    obj_4_el_.push_back(objects_[ele->dim()].size());
    objects_[ele->dim()].push_back(l);
  }
}


/**
 * The main functionality of the class creates duplicate nodes
 * at interfaces with fractures.
 * 
 * For each node:
 * 1) Create the lists of elements (node_elements) sharing the node.
 * 2) Divide node_elements into groups (components) connected by edges.
 * 3) If more than one component was created, create for each extra component a new node
 *    and update its index in the elements from this component.
 *    Also set node_dim_ (dimension of elements using the node).
 */
void DuplicateNodes::duplicate_nodes()
{
  START_TIMER("duplicate_nodes");
  for ( auto n : mesh_->node_range() ){
    // create list of adjacent elements from mesh_->node_elements[n.idx()]
    std::list<unsigned int> node_elements;
    std::copy( mesh_->node_elements()[n.idx()].begin(), mesh_->node_elements()[n.idx()].end(), std::back_inserter( node_elements ) );
    
    // divide node_elements into components connected by node_edges
    std::list<std::set<unsigned int> > components;
    while (node_elements.size() > 0)
    {
      // create component containing the first element in node_elements
      std::queue<unsigned int> q;
      q.push(*node_elements.begin());
      node_elements.erase(node_elements.begin());
      std::set<unsigned int> component;
      while (q.size() > 0) {
        auto elem = mesh_->element_accessor(q.front());
        component.insert(elem.idx());
        // add to queue adjacent elements sharing one of node_edges
        for (unsigned int sid=0; sid<elem->n_sides(); ++sid) {
          auto side = elem.side(sid);
          for (unsigned int esid=0; esid < side->edge().n_sides(); ++esid) {
            auto adj_el_idx = side->edge().side(esid)->elem_idx();
            if (adj_el_idx != elem.idx())
            {
              auto it = std::find(node_elements.begin(), node_elements.end(), adj_el_idx);
              if (it != node_elements.end()) {
                node_elements.erase(it);
                q.push(adj_el_idx);
              }
            }
          }
        }
        q.pop();
      }
      components.push_back(component);
    }
    
    // For the first component leave the current node and set its node_dim_.
    // If there are more components, then for each additional component
    // create a duplicate of the node and update the affected mesh objects.
    if (components.size() > 0) {
      node_dim_[n.idx()] = mesh_->element_accessor(*(components.begin()->begin())).dim();
      components.erase(components.begin());
    }
    for (auto component : components) {
      unsigned int new_nid = n_duplicated_nodes_++;
      node_dim_.push_back(mesh_->element_accessor(*component.begin()).dim());
      for (auto el_idx : component) {
        // After we have complete graph tetrahedron-triangles-lines-points,
        // the updating of duplicated nodes will have to change.
        MeshObject *elem = &objects_[mesh_->element_accessor(el_idx).dim()][obj_4_el_[el_idx]];
        for (unsigned int i=0; i<mesh_->element_accessor(el_idx).dim()+1; i++)
          if (elem->nodes[i] == n.idx())
            elem->nodes[i] = new_nid;
        for (unsigned int fi=0; fi<mesh_->element_accessor(el_idx).dim()+1; fi++)
          for (unsigned int ni=0; ni<mesh_->element_accessor(el_idx).dim(); ni++)
            if (elem->faces[fi]->nodes[ni] == n.idx())
              elem->faces[fi]->nodes[ni] = new_nid;
      }
    }
  }
  END_TIMER("duplicate_nodes");
}




