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
 * @file    mesh_tree.h
 * @brief   
 */

#ifndef MESH_TREE_H
#define MESH_TREE_H




class Mesh;

template<int dim>
class MeshObject {
public:
  
  MeshObject();
  ~MeshObject() {}
  
  MeshObject<dim-1> *faces[dim+1];
  unsigned int nodes[dim+1];
};


/**
 * Class MeshTree constructs the graph structure of elements, their faces and nodes
 * without any other data such as coordinates or region numbers. The nodes are then
 * duplicated where the elements are separated by fractures (elements of lower dim.).
 * 
 * The structure is used in DOF handler to distribute dofs for FE spaces sharing dofs
 * between elements.
 * 
 * TODO: Currently we do not create faces of faces, i.e. 1d edges of tetrahedra are not
 * available. This should be implemented if we need to distribute dofs on edges.
 * 
 */
class MeshTree {
public:
  
  MeshTree(Mesh *mesh);
  
  Mesh *mesh() const { return mesh_; }
  
  const std::vector<unsigned int> &nodes() const { return nodes_; }
  const std::vector<unsigned int> &node_dim() const { return node_dim_; }
  
  const std::vector<MeshObject<0> > &points() const { return points_; }
  const std::vector<MeshObject<1> > &lines() const { return lines_; }
  const std::vector<MeshObject<2> > &triangles() const { return triangles_; }
  const std::vector<MeshObject<3> > &tetras() const { return tetras_; }
  
  const std::vector<unsigned int> &obj_4_el() const { return obj_4_el_; }
  const std::vector<unsigned int> &obj_4_edg() const { return obj_4_edg_; }
  const std::vector<unsigned int> &obj_4_node() const { return obj_4_node_; }
  
  
private:
  
  void init_nodes();
  void init_from_edges();
  void init_from_elements();
  void duplicate_nodes();
  
  Mesh *mesh_;
  
  std::vector<unsigned int> nodes_;
  std::vector<unsigned int> node_dim_;
  
  std::vector<MeshObject<0> > points_;
  std::vector<MeshObject<1> > lines_;
  std::vector<MeshObject<2> > triangles_;
  std::vector<MeshObject<3> > tetras_;
  
  std::vector<unsigned int> obj_4_el_;
  std::vector<unsigned int> obj_4_edg_;
  std::vector<unsigned int> obj_4_node_;
};


#endif // MESH_TREE_H