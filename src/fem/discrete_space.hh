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
 * @file    discrete_space.hh
 * @brief   Declaration of class which provides the finite element for every mesh cell.
 * @author  Jan Stebel
 */

#ifndef DISCRETE_SPACE_HH_
#define DISCRETE_SPACE_HH_

#include "mesh/accessors.hh"


template<unsigned int dim> class FiniteElement;
class Mesh;


/**
 * Abstract class for definition of finite element functions on the mesh.
 * This should include
 * - simple FE spaces using the same finite element on
 *   all mesh elements with the same dimension,
 * - p-refined spaces using variable FE order,
 * - XFEM with arbitrary shape functions on every element.
 * 
 */
class DiscreteSpace {
public:
    
  /// Number of dofs associated to node. @p nid is the node index in the mesh tree.
  virtual unsigned int n_node_dofs(unsigned int nid) const = 0;

  /// Number of dofs associated to edge.
  virtual unsigned int n_edge_dofs(const Edge &edge) const = 0;
  
  /// Number of dofs associated to element (not shared by adjacent elements).
  virtual unsigned int n_elem_dofs(const ElementAccessor<3> &cell) const = 0;
  
  /// Number of dofs associated to generalized n-face (node, line, triangle or tetrahedron).
  template<unsigned int dim>
  unsigned int n_face_dofs(unsigned int face_id)
  {
    ASSERT(false).error("Not implemented.");
    return 0;
  }
  
  /// Return finite element object for given element.
  template<unsigned int dim>
  FiniteElement<dim> *fe(const ElementAccessor<3> &) const;

  /// Destructor.
  virtual ~DiscreteSpace() {};
    

protected:
  
  /// Constructor.
  DiscreteSpace(Mesh *mesh)
  : mesh_(mesh) {}
  
  virtual FiniteElement<0> *fe0d(const ElementAccessor<3> &) const = 0;
  virtual FiniteElement<1> *fe1d(const ElementAccessor<3> &) const = 0;
  virtual FiniteElement<2> *fe2d(const ElementAccessor<3> &) const = 0;
  virtual FiniteElement<3> *fe3d(const ElementAccessor<3> &) const = 0;
  
  Mesh *mesh_;

};



/**
 * Implementation of DiscreteSpace when all elements have the same FiniteElement.
 */
class EqualOrderDiscreteSpace : public DiscreteSpace {
public:
  
  EqualOrderDiscreteSpace(Mesh *mesh, FiniteElement<0> *fe0, FiniteElement<1> *fe1, FiniteElement<2> *fe2, FiniteElement<3> *fe3)
  : DiscreteSpace(mesh), fe0_(fe0), fe1_(fe1), fe2_(fe2), fe3_(fe3) {}
  
  unsigned int n_elem_dofs(const ElementAccessor<3> &cell) const override;
  
  unsigned int n_edge_dofs(const Edge &edge) const override;
  
  unsigned int n_node_dofs(unsigned int nid) const override;
  
  FiniteElement<0> *fe0d(const ElementAccessor<3> &) const override { return fe0_; }
  FiniteElement<1> *fe1d(const ElementAccessor<3> &) const override { return fe1_; }
  FiniteElement<2> *fe2d(const ElementAccessor<3> &) const override { return fe2_; }
  FiniteElement<3> *fe3d(const ElementAccessor<3> &) const override { return fe3_; }
  
  
private:
  
  FiniteElement<0> *fe0_;
  FiniteElement<1> *fe1_;
  FiniteElement<2> *fe2_;
  FiniteElement<3> *fe3_;
  
};







#endif /* DISCRETE_SPACE_HH_ */
