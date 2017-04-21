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
 * @file    fe_system.cc
 * @brief   Class FESystem for compound finite elements.
 * @author  Jan Stebel
 */

#include "fem/fe_system.hh"




template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::FESystem(std::vector<const FiniteElement<dim,spacedim> &> finite_elements)
  : fe_(finite_elements)
{
  this->init();

  this->is_primitive_ = false;
  
  unsigned int fe_index = 0;
  unsigned int comp_offset = 0;
  for (auto fe : fe_)
  {
    number_of_dofs += fe.n_dofs();
    this->n_components_ += fe.n_components();

    for (int i=0; i<=dim; ++i)
    {
      number_of_single_dofs[i] += fe.n_object_dofs(i, DOF_SINGLE);
      number_of_pairs[i] += fe.n_object_dofs(i, DOF_PAIR);
      number_of_triples[i] += fe.n_object_dofs(i, DOF_TRIPLE);
      number_of_sextuples[i] += fe.n_object_dofs(i, DOF_SEXTUPLE);
    }

    for (int i=0; i<fe.generalized_support_points().size(); i++)
        this->generalized_support_points.push_back(fe.get_generalized_support_points()[i]);

    if (fe.polynomial_order() > order) order = fe.polynomial_order();
    
    for (int i=0; i<fe.n_dofs(); ++i)
      fe_dof_indices_.push_back(DofComponentData(fe_index, i, comp_offset));
    
    fe_index++;
    comp_offset += fe.n_components();
  }
  
  if (this->is_primitive_)
  {
    double dof_index = 0;
    for (auto fe : fe_)
    {
      for (int i=0; i<fe.n_dofs(); ++i)
        this->component_indices_.push_back(fe_dof_indices_[dof_index++].fe_index);
    }
  } else {
    double dof_index = 0;
    comp_offset = 0;
    for (auto fe : fe_)
    {
      for (int i=0; i<fe.n_dofs(); ++i)
      {
        std::vector<bool> nonzeros(this->n_components_, false);
        std::copy_n(fe.get_nonzero_components(fe_dof_indices_[dof_index++].basis_index), fe.n_components(), &nonzeros[comp_offset]);
        this->nonzero_components_.push_back(nonzeros);
      }
      comp_offset += fe.n_components();
    }
  }

  this->compute_node_matrix();
}

template<unsigned int dim, unsigned int spacedim>
double FESystem<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  return fe_[fe_dof_indices_[i].first].basis_value(fe_dof_indices_[i].second, p);
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FESystem<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  return fe_[fe_dof_indices_[i].first].basis_grad(fe_dof_indices_[i].second, p);
}

template<unsigned int dim, unsigned int spacedim>
double FESystem<dim,spacedim>::basis_value_component(const unsigned int i, 
                                                     const arma::vec::fixed<dim> &p, 
                                                     const unsigned int comp) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  
  if (comp-fe_dof_indices_[i].component_offset >= 0 && 
      comp-fe_dof_indices_[i].component_offset < fe_[fe_dof_indices_[i].fe_index].n_components())
    return fe_[fe_dof_indices_[i].fe_index].basis_value_component(fe_dof_indices_[i].basis_index, p, comp-fe_dof_indices_[i].component_offset);
  
  return 0;
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FESystem<dim,spacedim>::basis_grad_component(const unsigned int i, 
                                                                   const arma::vec::fixed<dim> &p, 
                                                                   const unsigned int comp) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  
  if (comp-fe_dof_indices_[i].component_offset >= 0 && 
      comp-fe_dof_indices_[i].component_offset < fe_[fe_dof_indices_[i].fe_index].n_components())
    return fe_[fe_dof_indices_[i].fe_index].basis_grad_component(fe_dof_indices_[i].basis_index, p, comp-fe_dof_indices_[i].component_offset);
  
  return 0;
}

template<unsigned int dim, unsigned int spacedim>
void FESystem<dim,spacedim>::compute_node_matrix()
{
  OLD_ASSERT_EQUAL(this->get_generalized_support_points().size(), number_of_dofs);

  this->node_matrix.resize(number_of_dofs, number_of_dofs);

  unsigned int offset = 0;
  for (unsigned int i = 0; i < fe_.size(); i++)
  {
    this->node_matrix.submat(offset, offset, offset+fe_[i].n_dofs(), offset+fe_[i].n_dofs())
      = fe_[i].node_matrix;
      
    offset += fe_[i].n_dofs();
  }
}

template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::~FESystem()
{}










