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
#include "system/global_defs.h"
#include "quadrature/quadrature.hh"

using namespace std;



template<unsigned int dim,unsigned int spacedim>
std::vector<FiniteElement<dim,spacedim> > expand_fe_vector(const std::vector<std::pair<FiniteElement<dim,spacedim>, unsigned int> > &fe)
{
  std::vector<FiniteElement<dim,spacedim> > finite_elements;
  
  for (auto fe_pair : fe)
    for (unsigned int i=0; i<fe_pair.second; i++)
      finite_elements.push_back(fe_pair.first);

  return finite_elements;
}



template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::FESystem(FiniteElement<dim,spacedim> *fe, unsigned int n)
  : FiniteElement<dim,spacedim>()
{
  std::shared_ptr<FiniteElement<dim,spacedim> > fe_ptr(fe);
  fe_ = std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > >(n, fe_ptr);
  initialize();
}



template<unsigned int dim, unsigned int spacedim>
void FESystem<dim,spacedim>::initialize()
{
  this->is_primitive_ = false;
  n_components_ = 0;
  
  std::vector<unsigned int> node_basis[dim+1];
  std::vector<unsigned int> cell_basis;
  
  unsigned int fe_index = 0;
  unsigned int comp_offset = 0;
  unsigned basis_offset = 0;
  for (auto fe : fe_)
  {
    number_of_dofs += fe->n_dofs();
    n_components_ += fe->n_components();

    for (int i=0; i<=dim; ++i)
    {
      number_of_single_dofs[i] += fe->n_object_dofs(i, DOF_SINGLE);
      number_of_pairs[i] += fe->n_object_dofs(i, DOF_PAIR);
      number_of_triples[i] += fe->n_object_dofs(i, DOF_TRIPLE);
      number_of_sextuples[i] += fe->n_object_dofs(i, DOF_SEXTUPLE);
    }

    for (int i=0; i<fe->get_generalized_support_points().size(); i++)
        generalized_support_points.push_back(fe->get_generalized_support_points()[i]);

    if (fe->polynomial_order() > order) order = fe->polynomial_order();
    
    for (int i=0; i<fe->n_dofs(); ++i)
      fe_dof_indices_.push_back(DofComponentData(fe_index, i, comp_offset));
    
    for (unsigned int i=0; i<fe->n_object_dofs(0, DOF_SINGLE); i++)
      node_basis[i % (dim+1)].push_back(basis_offset+i);
    
    unsigned int n_cell_dofs = fe->n_object_dofs(dim, DOF_SINGLE)
                              +fe->n_object_dofs(dim, DOF_PAIR)
                              +fe->n_object_dofs(dim, DOF_TRIPLE)
                              +fe->n_object_dofs(dim, DOF_SEXTUPLE);
    for (unsigned int i=0; i<n_cell_dofs; i++)
      cell_basis.push_back(basis_offset+fe->n_object_dofs(0,DOF_SINGLE)+i);
    
    fe_index++;
    comp_offset += fe->n_components();
    basis_offset += fe->n_dofs();
  }
  
  for (unsigned int i=0; i<dim+1; i++)
    dof_basis_.insert(dof_basis_.end(), node_basis[i].begin(), node_basis[i].end());
  dof_basis_.insert(dof_basis_.end(), cell_basis.begin(), cell_basis.end());
  
  if (this->is_primitive_)
  {
    double dof_index = 0;
    for (auto fe : fe_)
    {
      for (int i=0; i<fe->n_dofs(); ++i)
        this->component_indices_.push_back(fe_dof_indices_[dof_index++].fe_index);
    }
  } else {
    double dof_index = 0;
    comp_offset = 0;
    for (auto fe : fe_)
    {
      for (int i=0; i<fe->n_dofs(); ++i)
      {
        std::vector<bool> nonzeros(n_components_, false);
        for (unsigned int c=0; c<fe->n_components(); c++)\
          nonzeros[comp_offset+c] = fe->get_nonzero_components(fe_dof_indices_[dof_index++].basis_index)[c];
        this->nonzero_components_.push_back(nonzeros);
      }
      comp_offset += fe->n_components();
    }
  }

  compute_node_matrix();
}



template<unsigned int dim, unsigned int spacedim>
double FESystem<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  return fe_[fe_dof_indices_[dof_basis_[i]].fe_index]->basis_value(fe_dof_indices_[dof_basis_[i]].basis_index, p);
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FESystem<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  return fe_[fe_dof_indices_[dof_basis_[i]].fe_index]->basis_grad(fe_dof_indices_[dof_basis_[i]].basis_index, p);
}

template<unsigned int dim, unsigned int spacedim>
double FESystem<dim,spacedim>::basis_value_component(const unsigned int i, 
                                                     const arma::vec::fixed<dim> &p, 
                                                     const unsigned int comp) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  
  unsigned int bi = dof_basis_[i];
  int l_comp = comp-fe_dof_indices_[bi].component_offset;
  if (l_comp >= 0 && l_comp < fe_[fe_dof_indices_[bi].fe_index]->n_components())
    if (fe_[fe_dof_indices_[bi].fe_index]->n_components() == 1)
      return fe_[fe_dof_indices_[bi].fe_index]->basis_value(fe_dof_indices_[bi].basis_index, p);
    else
      return fe_[fe_dof_indices_[bi].fe_index]->basis_value_component(fe_dof_indices_[bi].basis_index, p, l_comp);
  
  return 0;
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FESystem<dim,spacedim>::basis_grad_component(const unsigned int i, 
                                                                   const arma::vec::fixed<dim> &p, 
                                                                   const unsigned int comp) const
{
  OLD_ASSERT(i <= number_of_dofs, "Index of basis function is out of range.");
  
  unsigned int bi = dof_basis_[i];
  int l_comp = comp-fe_dof_indices_[bi].component_offset;
  if (l_comp >= 0 && l_comp < fe_[fe_dof_indices_[bi].fe_index]->n_components())
    if (fe_[fe_dof_indices_[bi].fe_index]->n_components() == 1)
      return fe_[fe_dof_indices_[bi].fe_index]->basis_grad(fe_dof_indices_[bi].basis_index, p);
    else
      return fe_[fe_dof_indices_[bi].fe_index]->basis_grad_component(fe_dof_indices_[bi].basis_index, p, l_comp);
  
//   return 0;
}


template<unsigned int dim, unsigned int spacedim> inline
UpdateFlags FESystem<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    for (auto fe : fe_)
      f |= fe->update_each(flags);

    return f;
}


template<unsigned int dim, unsigned int spacedim>
void FESystem<dim,spacedim>::compute_node_matrix()
{
  OLD_ASSERT_EQUAL(this->get_generalized_support_points().size(), number_of_dofs);

  this->node_matrix.resize(number_of_dofs, number_of_dofs);

  unsigned int offset = 0;
  for (unsigned int i = 0; i < fe_.size(); i++)
  {
    this->node_matrix.submat(offset, offset, offset+fe_[i]->n_dofs()-1, offset+fe_[i]->n_dofs()-1)
      = fe_[i]->node_matrix;
      
    offset += fe_[i]->n_dofs();
  }
}


template<unsigned int dim, unsigned int spacedim>
FEInternalData *FESystem<dim,spacedim>::initialize(const Quadrature<dim> &q, UpdateFlags flags)
{
    FEInternalData *data = new FEInternalData;
    std::vector<FEInternalData *> fe_data;
    
    if ((flags & update_values) || (flags & update_gradients))
      for (auto fe : fe_)
        fe_data.push_back(fe->initialize(q, flags));

    if (flags & update_values)
    {
        data->basis_vectors.resize(q.size(), std::vector<arma::vec>(number_of_dofs, arma::vec(n_components_)));
        
        unsigned int comp_offset = 0;
        unsigned int dof_offset = 0;
        for (unsigned int f=0; f<fe_.size(); f++)
        {
          if (fe_[f]->n_components() == 1)
          {
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                data->basis_vectors[i][dof_offset+n][comp_offset] = fe_data[f]->basis_values[i][n];
          }
          else
          {
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                data->basis_vectors[i][dof_offset+n].subvec(comp_offset,comp_offset+fe_[f]->n_components()-1) = fe_data[f]->basis_vectors[i][n];
          }
          comp_offset += fe_[f]->n_components();
          dof_offset += fe_[f]->n_dofs();
        }
    }

    if (flags & update_gradients)
    {
        data->basis_grad_vectors.resize(q.size(), std::vector<arma::mat>(number_of_dofs, arma::mat(n_components_,dim)));
        
        unsigned int comp_offset = 0;
        unsigned int dof_offset = 0;
        for (unsigned int f=0; f<fe_.size(); f++)
        {
          if (fe_[f]->n_components() == 1)
          {
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                data->basis_grad_vectors[i][dof_offset+n].row(comp_offset) = fe_data[f]->basis_grads[i].row(n);
          }
          else
          {
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                data->basis_grad_vectors[i][dof_offset+n].submat(comp_offset,0,comp_offset+fe_[f]->n_components()-1,dim-1) = fe_data[f]->basis_grad_vectors[i][n];
          }
          comp_offset += fe_[f]->n_components();
          dof_offset += fe_[f]->n_dofs();
        }
    }
    
    for (auto d : fe_data) delete d;

    return data;
}

template<unsigned int dim, unsigned int spacedim> inline
void FESystem<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &q,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    if ((fv_data.update_flags & update_values) || (fv_data.update_flags & update_gradients))
    {
      std::vector<std::vector<double> > shape_values(q.size(), std::vector<double>(fv_data.shape_values[0].size()));
      std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_gradients(q.size(), std::vector<arma::vec::fixed<spacedim> >(fv_data.shape_values[0].size()));
      unsigned int comp_offset = 0;
      unsigned int dof_offset = 0;
      unsigned int shape_offset = 0;
      for (unsigned int f=0; f<fe_.size(); f++)
      {
        FEInternalData d;
        FEValuesData<dim,spacedim> sub_fv_data;
        
        if (fv_data.update_flags & update_values)
        {
          if (fe_[f]->n_components() == 1)
          {
            d.basis_values.resize(q.size(), arma::vec(fe_[f]->n_dofs()));
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                d.basis_values[i].row(n) = data.basis_vectors[i][dof_offset+n].row(comp_offset);
          }
          else
          {
            d.basis_vectors.resize(q.size(), std::vector<arma::vec>(fe_[f]->n_dofs(), arma::vec(fe_[f]->n_components())));
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                d.basis_vectors[i][n] = data.basis_vectors[i][dof_offset+n].subvec(comp_offset,comp_offset+fe_[f]->n_components()-1);
          }
        }
        
        if (fv_data.update_flags & update_gradients)
        {
          if (fe_[f]->n_components() == 1)
          {
            d.basis_grads.resize(q.size(), arma::mat(fe_[f]->n_dofs(),dim));
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                d.basis_grads[i].row(n) = data.basis_grad_vectors[i][dof_offset+n].row(comp_offset);
          }
          else
          {
            d.basis_grad_vectors.resize(q.size(), std::vector<arma::mat>(fe_[f]->n_dofs(), arma::mat(fe_[f]->n_components(),dim)));
            for (unsigned int i=0; i<q.size(); i++)
              for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
                d.basis_grad_vectors[i][n] = data.basis_grad_vectors[i][dof_offset+n].submat(comp_offset,0,comp_offset+fe_[f]->n_components()-1,dim-1);
          }
        }
        
        sub_fv_data.allocate(q.size(), fv_data.update_flags, fe_[f]->n_dofs()*fe_[f]->n_components());
        sub_fv_data.jacobians = fv_data.jacobians;
        sub_fv_data.inverse_jacobians = fv_data.inverse_jacobians;
        sub_fv_data.determinants = fv_data.determinants;
        fe_[f]->fill_fe_values(q, d, sub_fv_data);
        
        if (fv_data.update_flags & update_values)
        {
          for (unsigned int i=0; i<q.size(); i++)
            for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
              for (unsigned int c=0; c<fe_[f]->n_components(); c++)
                shape_values[i][shape_offset+n_components_*n+comp_offset+c] = sub_fv_data.shape_values[i][n];
        }
        
        if (fv_data.update_flags & update_gradients)
        {
          for (unsigned int i=0; i<q.size(); i++)
            for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
              for (unsigned int c=0; c<fe_[f]->n_components(); c++)
                shape_gradients[i][shape_offset+n_components_*n+comp_offset+c] = sub_fv_data.shape_gradients[i][n];
        }
        
        dof_offset += fe_[f]->n_dofs();
        comp_offset += fe_[f]->n_components();
        shape_offset += fe_[f]->n_dofs()*n_components_;
      }
    
      // reorder values for compatibility with dof handler
      for (unsigned int i=0; i<q.size(); i++)
      {
        for (unsigned int n=0; n<number_of_dofs; n++)
          for (unsigned int c=0; c<n_components_; c++)
          {
            if (fv_data.update_flags & update_values)
              fv_data.shape_values[i][n*n_components_+c] = shape_values[i][dof_basis_[n]*n_components_+c];
            if (fv_data.update_flags & update_gradients)
              fv_data.shape_gradients[i][n*n_components_+c] = shape_gradients[i][dof_basis_[n]*n_components_+c];
          }
      }
      
    }
}


template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::~FESystem()
{}






template class FESystem<1,3>;
template class FESystem<2,3>;
template class FESystem<3,3>;



