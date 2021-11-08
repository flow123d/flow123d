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
 * @file    duplicate_nodes.h
 * @brief   
 */

#ifndef DUPLICATE_NODES_H
#define DUPLICATE_NODES_H




class MeshBase;

/**
 * Class representing an n-face (node, line, triangle or tetrahedron) in the mesh.
 * Currently, the object stores information about its 1-codimensional faces
 * and nodes.
 * 
 * In the future we should think about complete graph, where every n-dimensional object
 * has information only about its (n-1)-dimensional faces.
 */
class MeshObject {
public:
  
  MeshObject(unsigned int dim);
  ~MeshObject() {}
  
  /// Indices of (dim-1)-dimensional faces within the global vector of objects with given dimension.
  unsigned int faces[4];
  
  /// Indices of nodes.
  unsigned int nodes[4];
  
  /// Dimension of n-face.
  unsigned int dim_;
};


/**
 * Class DuplicateNodes constructs the graph structure of elements, their faces and nodes
 * without any other data such as coordinates or region numbers. The nodes are then
 * duplicated where the elements are separated by fractures (elements of lower dim.).
 * E.g.:
 * 
 * Consider a domain containing fracture like this:
 * +-----+-----+
 * |     |     |
 * |     +     |
 * |           |
 * |           |
 * |           |
 * +-----------+
 * 
 * with the mesh nodes numbered as follows:
 * 0-----1-----2
 * |     |     |
 * |     3     |
 * |           |
 * |           |
 * |           |
 * 4-----------5
 * 
 * The class duplicates node 1 and node 3 because in these nodes FE functions can be discontinuous
 * (schematic view):
 * 0-6   1   7-2
 * |  \  |  /  |
 * |   \ 3 /   |
 * |    \ /    |
 * |     8     |
 * |           |
 * 4-----------5
 * 
 * Now, nodes 1,6,7 have the same coordinates but belong to different groups of elements
 * and analogously nodes 3 and 8.
 * 
 * 
 * The structure is used in DOF handler to distribute dofs for FE spaces sharing dofs
 * between elements.
 * 
 * TODO: Currently we do not create faces of faces, i.e. 1d edges of tetrahedra are not
 * available. This should be implemented if we need to distribute dofs on edges.
 * 
 */
class DuplicateNodes {
public:
  
  DuplicateNodes(MeshBase *mesh);

  // Getters (see private section for explanation).  
  MeshBase *mesh() const { return mesh_; }
  
  unsigned int n_nodes() const { return n_duplicated_nodes_; }
  const std::vector<unsigned int> &node_dim() const { return node_dim_; }
  
  const std::vector<MeshObject> &objects(unsigned int dim) const
  { return objects_[dim]; }
  
  const std::vector<unsigned int> &obj_4_el() const { return obj_4_el_; }
  const std::vector<unsigned int> &obj_4_edg() const { return obj_4_edg_; }
  
  
private:
  
  /// Initialize the vector of nodes from mesh.
  void init_nodes();
  
  /// Initialize objects from mesh edges.
  void init_from_edges();
  
  /// Initialize objects from mesh elements.
  void init_from_elements();
  
  /// Duplicate nodes that are lying on interfaces with fractures.
  void duplicate_nodes();
  
  
  /// The mesh object.
  MeshBase *mesh_;
  
  /// Number of nodes (including duplicated ones).
  unsigned int n_duplicated_nodes_;
  
  /// Vector of space dimensions of elements using the particular duplicated node.
  std::vector<unsigned int> node_dim_;
  
  /// Array of n-faces by their dimension.
  std::vector<MeshObject> objects_[4];
  
  /** Vector of object indices for each mesh element.
   * For an element with index el_idx, obj_4_el_[el_idx] is the index of the corresponding object
   * in the vector tetras_/triangles_/lines_/points_, depending on the element dimension.
   */
  std::vector<unsigned int> obj_4_el_;
  
  /// Vector of object indices for each mesh edge.
  std::vector<unsigned int> obj_4_edg_;
};


#endif // DUPLICATE_NODES_H
