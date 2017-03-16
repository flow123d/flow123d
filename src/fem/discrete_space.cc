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
 * @file    dofhandler.hh
 * @brief   Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 * @author  Jan Stebel
 */

#include "fem/discrete_space.hh"
#include "mesh/mesh.h"
#include "mesh/mesh_tree.h"



unsigned int EqualOrderDiscreteSpace::n_elem_dofs(const ElementFullIter &cell) const
{
  switch (cell->dim())
  {
    case 1: return fe1_->n_object_dofs(cell->dim(), DofMultiplicity::DOF_SINGLE);
    case 2: return fe2_->n_object_dofs(cell->dim(), DofMultiplicity::DOF_SINGLE);
    case 3: return fe3_->n_object_dofs(cell->dim(), DofMultiplicity::DOF_SINGLE);
  }
}


unsigned int EqualOrderDiscreteSpace::n_node_dofs(unsigned int nid) const
{
  unsigned int dim = mesh_->tree->node_dim()[nid];
  switch (dim)
  {
    case 1: return fe1_->n_object_dofs(0, DofMultiplicity::DOF_SINGLE)/(dim+1);
    case 2: return fe2_->n_object_dofs(0, DofMultiplicity::DOF_SINGLE)/(dim+1);
    case 3: return fe3_->n_object_dofs(0, DofMultiplicity::DOF_SINGLE)/(dim+1);
  }
}



template<> FiniteElement<1,3> *DiscreteSpace::fe(const ElementFullIter &cell) const { return fe1d(cell); }
template<> FiniteElement<2,3> *DiscreteSpace::fe(const ElementFullIter &cell) const { return fe2d(cell); }
template<> FiniteElement<3,3> *DiscreteSpace::fe(const ElementFullIter &cell) const { return fe3d(cell); }





