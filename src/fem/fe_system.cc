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



template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::FESystem(std::shared_ptr<FiniteElement<dim,spacedim> > fe, FEType t)
{
  OLD_ASSERT(fe->is_primitive(), "FE vector or tensor can olny by created from primitive FE.");
  OLD_ASSERT(t == FEType::FEVector || t == FEType::FETensor, "This constructor can be used only for vectors or tensors.");
  
  FiniteElement<dim,spacedim>::init(false, t);
  
  if (t == FEType::FEVector)
    fe_ = std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > >(dim, fe);
  else
    fe_ = std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > >(dim*dim, fe);
  
  initialize();
}


template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::FESystem(const std::shared_ptr<FiniteElement<dim,spacedim> > &fe, unsigned int n)
{
  FiniteElement<dim,spacedim>::init(false, FEMixedSystem);
  fe_ = std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > >(n, fe);
  initialize();
}


template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::FESystem(std::vector<std::shared_ptr<FiniteElement<dim,spacedim> > > fe)
{
  FiniteElement<dim,spacedim>::init(false, FEMixedSystem);
  for (std::shared_ptr<FiniteElement<dim,spacedim> > fe_object : fe)
    fe_.push_back(fe_object);
  initialize();
}



template<unsigned int dim, unsigned int spacedim>
void FESystem<dim,spacedim>::initialize()
{
  unsigned int fe_index = 0;
  unsigned int comp_offset = 0;
  n_components_ = 0;
  // for each base FE add components, support points, and other 
  // information to the system
  for (auto fe : fe_)
  {
    n_components_ += fe->n_components();
    
    switch (fe->type_)
    {
      case FEType::FEScalar:
        scalar_components_.push_back(comp_offset);
        break;
      case FEType::FEVector:
        vector_components_.push_back(comp_offset);
        break;
      default:
        OLD_ASSERT(false, "Not implemented.");
        break;
    }

    for (int i=0; i<fe->n_dofs(); ++i)
      fe_dof_indices_.push_back(DofComponentData(fe_index, i, comp_offset));
    
    fe_index++;
    comp_offset += fe->n_components();
  }
  
  double dof_index = 0;
  comp_offset = 0;
  for (auto fe : fe_)
  {
      for (unsigned int i=0; i<fe->n_dofs(); i++)
      {
          arma::vec coefs(n_components_);
          coefs.subvec(comp_offset, comp_offset+fe->dof(i).coefs.size()-1) = fe->dof(i).coefs;
          this->dofs_.push_back(Dof(fe->dof(i).dim, fe->dof(i).n_face_idx, fe->dof(i).coords, coefs, fe->dof(i).type));
          dof_index++;
      }
      comp_offset += fe->n_components();
  }
  
//   if (this->is_primitive_)
//   {
//     double dof_index = 0;
//     // add index of nonzero component for each dof in FESystem
//     for (auto fe : fe_)
//     {
//       for (int i=0; i<fe->n_dofs(); ++i)
//         this->component_indices_.push_back(fe_dof_indices_[dof_index++].fe_index);
//     }
//   } else {
    dof_index = 0;
    comp_offset = 0;
    // add footprint of nonzero components for each dof in FESystem
    for (auto fe : fe_)
    {
      for (int i=0; i<fe->n_dofs(); ++i)
      {
        std::vector<bool> nonzeros(n_components_, false);
        for (unsigned int c=0; c<fe->n_components(); c++)\
          nonzeros[comp_offset+c] = fe->get_nonzero_components(fe_dof_indices_[dof_index].basis_index)[c];
        this->nonzero_components_.push_back(nonzeros);
        dof_index++;
      }
      comp_offset += fe->n_components();
    }
//   }

  compute_node_matrix();
}



template<unsigned int dim, unsigned int spacedim>
double FESystem<dim,spacedim>::basis_value(const unsigned int i, 
                                           const arma::vec::fixed<dim> &p, 
                                           const unsigned int comp) const
{
  OLD_ASSERT(i <= this->dofs_.size(), "Index of basis function is out of range.");
  
  // component index in the base FE
  int l_comp = comp-fe_dof_indices_[i].component_offset;
  OLD_ASSERT(l_comp >= 0 && l_comp < fe_[fe_dof_indices_[i].fe_index]->n_components(),
    "Index of component is out of range.");

  return fe_[fe_dof_indices_[i].fe_index]->basis_value(fe_dof_indices_[i].basis_index, p, l_comp);
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FESystem<dim,spacedim>::basis_grad(const unsigned int i, 
                                                         const arma::vec::fixed<dim> &p, 
                                                         const unsigned int comp) const
{
  OLD_ASSERT(i <= this->dofs_.size(), "Index of basis function is out of range.");
  
  // component index in the base FE
  int l_comp = comp-fe_dof_indices_[i].component_offset;
  OLD_ASSERT(l_comp >= 0 && l_comp < fe_[fe_dof_indices_[i].fe_index]->n_components(),
    "Index of component is out of range.");
  
  return fe_[fe_dof_indices_[i].fe_index]->basis_grad(fe_dof_indices_[i].basis_index, p, l_comp);
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
  // form the node_matrix of the FESystem as block diagonal matrix
  // composed of node_matrices of each base FE class
  
  this->node_matrix.resize(this->dofs_.size(), this->dofs_.size());

  unsigned int offset = 0;
  for (unsigned int i = 0; i < fe_.size(); i++)
  {
    this->node_matrix.submat(offset, offset, offset+fe_[i]->n_dofs()-1, offset+fe_[i]->n_dofs()-1)
      = fe_[i]->node_matrix;
      
    offset += fe_[i]->n_dofs();
  }
}


template<unsigned int dim, unsigned int spacedim>
FEInternalData *FESystem<dim,spacedim>::initialize(const Quadrature<dim> &q)
{
  FEInternalData *data = new FEInternalData;
  std::vector<FEInternalData *> fe_data;
  
  // first initialize the base FE
  for (auto fe : fe_)
    fe_data.push_back(fe->initialize(q));

  // fill values of basis functions
  data->basis_vectors.resize(q.size(), std::vector<arma::vec>(this->dofs_.size(), arma::vec(n_components_)));
  
  unsigned int comp_offset = 0;
  unsigned int dof_offset = 0;
  for (unsigned int f=0; f<fe_.size(); f++)
  {
    if (fe_[f]->n_components() == 1)
    {
      // for scalar base FE copy only one value per point and dof
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          data->basis_vectors[i][dof_offset+n][comp_offset] = fe_data[f]->basis_values[i][n];
    }
    else
    {
      // for vector-valued base FE copy the values to subvector
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          data->basis_vectors[i][dof_offset+n].subvec(comp_offset,comp_offset+fe_[f]->n_components()-1) = fe_data[f]->basis_vectors[i][n];
    }
    comp_offset += fe_[f]->n_components();
    dof_offset += fe_[f]->n_dofs();
  }

  // fill gradients of basis functions
  data->basis_grad_vectors.resize(q.size(), std::vector<arma::mat>(this->dofs_.size(), arma::mat(n_components_,dim)));
  
  comp_offset = 0;
  dof_offset = 0;
  for (unsigned int f=0; f<fe_.size(); f++)
  {
    if (fe_[f]->n_components() == 1)
    {
      // for scalar base FE copy values to one row per point and dof
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          data->basis_grad_vectors[i][dof_offset+n].row(comp_offset) = fe_data[f]->basis_grads[i].row(n);
    }
    else
    {
      // for vector-valued base FE copy the values to submatrix
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          data->basis_grad_vectors[i][dof_offset+n].submat(comp_offset,0,comp_offset+fe_[f]->n_components()-1,dim-1) = fe_data[f]->basis_grad_vectors[i][n];
    }
    comp_offset += fe_[f]->n_components();
    dof_offset += fe_[f]->n_dofs();
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
  // This method essentially calls fill_fe_values for the FE classes of components
  // of the FESystem and copies the values to the appropriate structures for FESystem.
  if (!(fv_data.update_flags & update_values) && !(fv_data.update_flags & update_gradients)) return;
  
  std::vector<std::vector<double> > shape_values(q.size(), std::vector<double>(fv_data.shape_values[0].size(), 0.));
  std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_gradients(q.size(), std::vector<arma::vec::fixed<spacedim> >(fv_data.shape_values[0].size()));
  for (unsigned int k=0; k<q.size(); k++)
    for (unsigned int i=0; i<fv_data.shape_values[0].size(); i++)
      shape_gradients[k][i].zeros();
  
  unsigned int comp_offset = 0;
  unsigned int dof_offset = 0;
  unsigned int shape_offset = 0;
  for (unsigned int f=0; f<fe_.size(); f++)
  {
    // create new FEInternalData object for the base FE
    // and fill its values/gradients
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
    
    // fill fe_values for base FE
    sub_fv_data.allocate(q.size(), fv_data.update_flags, fe_[f]->n_dofs()*fe_[f]->n_components());
    sub_fv_data.jacobians = fv_data.jacobians;
    sub_fv_data.inverse_jacobians = fv_data.inverse_jacobians;
    sub_fv_data.determinants = fv_data.determinants;
    fe_[f]->fill_fe_values(q, d, sub_fv_data);
    
    // gather fe_values in vectors for FESystem
    if (fv_data.update_flags & update_values)
    {
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          for (unsigned int c=0; c<fe_[f]->n_components(); c++)
            fv_data.shape_values[i][shape_offset+n_components_*n+comp_offset+c] = sub_fv_data.shape_values[i][n*fe_[f]->n_components()+c];
    }
    
    if (fv_data.update_flags & update_gradients)
    {
      for (unsigned int i=0; i<q.size(); i++)
        for (unsigned int n=0; n<fe_[f]->n_dofs(); n++)
          for (unsigned int c=0; c<fe_[f]->n_components(); c++)
            fv_data.shape_gradients[i][shape_offset+n_components_*n+comp_offset+c] = sub_fv_data.shape_gradients[i][n*fe_[f]->n_components()+c];
    }
    
    dof_offset += fe_[f]->n_dofs();
    comp_offset += fe_[f]->n_components();
    shape_offset += fe_[f]->n_dofs()*n_components_;
  }
}


template<unsigned int dim, unsigned int spacedim>
FESystem<dim,spacedim>::~FESystem()
{}






template class FESystem<1,3>;
template class FESystem<2,3>;
template class FESystem<3,3>;



