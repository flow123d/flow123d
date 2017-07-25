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

#include "mesh/duplicate_nodes.h"
#include "mesh/mesh.h"




template<int dim>
MeshObject<dim>::MeshObject()
{
  for (unsigned int i=0; i<dim+1; i++)
    faces[i] = nullptr;
}


DuplicateNodes::DuplicateNodes(Mesh* mesh)
: mesh_(mesh)
{
  init_nodes();
  init_from_edges();
  init_from_elements();
  duplicate_nodes();
}



void DuplicateNodes::init_nodes()
{
  // initialize the nodes
  unsigned int nid=0;
  FOR_NODES(mesh_, n) {
    obj_4_node_.push_back(nid);
    nodes_.push_back(nid++);
  }
  node_dim_.resize(nodes_.size());
}



void DuplicateNodes::init_from_edges()
{
  // initialize the objects from edges
  FOR_EDGES(mesh_, edge) {
    switch (edge->side(0)->dim()) {
      case 0:
        {
          MeshObject<0> p;
          p.nodes[0] = obj_4_node_[mesh_->node_vector.index(edge->side(0)->node(0))];
          obj_4_edg_.push_back(points_.size());
          points_.push_back(p);
        }
        break;
      case 1:
        {
          MeshObject<1> l;
          for (unsigned int i=0; i<2; i++)
            l.nodes[i] = obj_4_node_[mesh_->node_vector.index(edge->side(0)->node(i))];
          obj_4_edg_.push_back(lines_.size());
          lines_.push_back(l);
        }
        break;
      case 2:
        {
          MeshObject<2> t;
          for (unsigned int i=0; i<3; i++)
            t.nodes[i] = obj_4_node_[mesh_->node_vector.index(edge->side(0)->node(i))];
          obj_4_edg_.push_back(triangles_.size());
          triangles_.push_back(t);
        }
        break;
    }
  }
}



void DuplicateNodes::init_from_elements()
{
  // initialize the objects from elements
  FOR_ELEMENTS(mesh_, ele) {
    switch (ele->dim()) {
      case 1:
        {
          MeshObject<1> l;
          for (unsigned int i=0; i<2; i++)
            l.nodes[i] = obj_4_node_[mesh_->node_vector.index(ele->node[i])];
          for (unsigned int i=0; i<2; i++)
            l.faces[i] = &points_[obj_4_edg_[ele->side(i)->edge_idx()]];
          obj_4_el_.push_back(lines_.size());
          lines_.push_back(l);
        }
        break;
      case 2:
        {
          MeshObject<2> tr;
          for (unsigned int i=0; i<3; i++)
            tr.nodes[i] = obj_4_node_[mesh_->node_vector.index(ele->node[i])];
          for (unsigned int i=0; i<3; i++)
            tr.faces[i] = &lines_[obj_4_edg_[ele->side(i)->edge_idx()]];
          obj_4_el_.push_back(triangles_.size());
          triangles_.push_back(tr);
        }
        break;
      case 3:
        {
          MeshObject<3> te;
          for (unsigned int i=0; i<4; i++)
            te.nodes[i] = obj_4_node_[mesh_->node_vector.index(ele->node[i])];
          for (unsigned int i=0; i<4; i++)
            te.faces[i] = &triangles_[obj_4_edg_[ele->side(i)->edge_idx()]];
          obj_4_el_.push_back(tetras_.size());
          tetras_.push_back(te);
        }
        break;
    }
  }
}


/**
 * The main functionality of the class creates duplicate nodes
 * at interfaces with fractures.
 * 
 * For each node:
 * 1) Create the lists of elements (node_elements) and edges (node_edges) sharing the node.
 * 2) Divide node_elements into groups (components) which share one of node_edges.
 * 3) If more than one component was created, create for each extra component a new node
 *    and update its index in the elements from this component.
 *    Also set node_dim_ (dimension of elements using the node).
 */
void DuplicateNodes::duplicate_nodes()
{
  FOR_NODES(mesh_, n){
    // create list of elements
    std::vector<unsigned int> node_elements = mesh_->node_elements[n.index()];
    
    // create list of edges sharing the node
    std::set<unsigned int> node_edges;
    for (unsigned int el_idx : node_elements)
    {
      FOR_ELEMENT_SIDES(&mesh_->element[el_idx], sid)
        node_edges.insert(mesh_->element[el_idx].side(sid)->edge_idx());
    }
    
    // divide node_elements into components connected by node_edges
    std::vector<std::set<unsigned int> > components;
    while (node_elements.size() > 0)
    {
      // create component containing the first element in node_elements
      std::queue<unsigned int> q;
      q.push(node_elements[0]);
      node_elements.erase(node_elements.begin());
      std::set<unsigned int> component;
      while (q.size() > 0) {
        auto elem = mesh_->element.full_iter(&mesh_->element[q.front()]);
        component.insert(elem.index());
        // add to queue adjacent elements sharing one of node_edges
        FOR_ELEMENT_SIDES(elem, sid) {
          auto side = elem->side(sid);
          if (node_edges.find(side->edge_idx()) != node_edges.end())
          {
            unsigned int esid;
            FOR_EDGE_SIDES(side->edge(), esid) {
              auto adj_el_idx = side->edge()->side(esid)->element().index();
              if (adj_el_idx != elem.index())
              {
                auto it = std::find(node_elements.begin(), node_elements.end(), adj_el_idx);
                if (it != node_elements.end()) {
                  node_elements.erase(it);
                  q.push(adj_el_idx);
                }
              }
            }
          }
        }
        q.pop();
      }
      components.push_back(component);
    }
    
    // for each component except the first one
    // create a duplicate of the node
    if (components.size() > 0) {
      node_dim_[obj_4_node_[n.index()]] = mesh_->element[*components[0].begin()].dim();
      components.erase(components.begin());
    }
    for (auto component : components) {
      unsigned int new_nid = nodes_.size();
      nodes_.push_back(new_nid);
      node_dim_.push_back(mesh_->element[*component.begin()].dim());
      for (auto el_idx : component) {
        switch (mesh_->element[el_idx].dim()) {
          case 1:
            {
              MeshObject<1> *elem = &lines_[obj_4_el_[el_idx]];
              for (unsigned int i=0; i<2; i++)
                if (elem->nodes[i] == n.index())
                  elem->nodes[i] = new_nid;
              for (unsigned int fi=0; fi<2; fi++)
                if (elem->faces[fi]->nodes[0] == n.index())
                  elem->faces[fi]->nodes[0] = new_nid;
            }
            break;
          case 2:
            {
              MeshObject<2> *elem = &triangles_[obj_4_el_[el_idx]];
              for (unsigned int i=0; i<3; i++)
                if (elem->nodes[i] == n.index())
                  elem->nodes[i] = new_nid;
              for (unsigned int fi=0; fi<3; fi++)
                for (unsigned int ni=0; ni<2; ni++)
                  if (elem->faces[fi]->nodes[ni] == n.index())
                    elem->faces[fi]->nodes[ni] = new_nid;
            }
            break;
          case 3:
            {
              MeshObject<3> *elem = &tetras_[obj_4_el_[el_idx]];
              for (unsigned int i=0; i<4; i++)
                if (elem->nodes[i] == n.index())
                  elem->nodes[i] = new_nid;
              for (unsigned int fi=0; fi<4; fi++)
                for (unsigned int ni=0; ni<3; ni++)
                  if (elem->faces[fi]->nodes[ni] == n.index())
                    elem->faces[fi]->nodes[ni] = new_nid;
            }
            break;
        }
      }
    }
  }
}




