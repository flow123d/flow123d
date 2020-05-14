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

#include "mesh/accessors.hh"
#include "fem/fe_system.hh"
#include "system/global_defs.h"
#include "quadrature/quadrature.hh"

using namespace std;


unsigned int count_components(const vector<shared_ptr<FunctionSpace> > &fs_vector)
{
    unsigned int n_comp = 0;
    for (auto fs : fs_vector)
        n_comp += fs->n_components();
    
    return n_comp;
}


unsigned int check_spacedim(const vector<shared_ptr<FunctionSpace> > &fs_vector)
{
    ASSERT_DBG(fs_vector.size() > 0);
    unsigned int space_dim = fs_vector[0]->space_dim();
    for (auto fs : fs_vector)
        ASSERT_DBG(fs->space_dim() == space_dim).error("FunctionSpace space_dim mismatch.");
    
    return space_dim;
}


FESystemFunctionSpace::FESystemFunctionSpace(const vector<shared_ptr<FunctionSpace> > &fs_vector)
    : FunctionSpace(check_spacedim(fs_vector), count_components(fs_vector)),
      fs_(fs_vector)
{
    dim_ = 0;
    unsigned int fe_index = 0;
    unsigned int comp_offset = 0;
    for (auto fs : fs_vector)
    {
        for (unsigned int i=0; i<fs->dim(); i++)
            dof_indices_.push_back(DofComponentData(fe_index, i, comp_offset));
        
        dim_ += fs->dim();
        fe_index++;
        comp_offset += fs->n_components();
    }
}


double FESystemFunctionSpace::basis_value(unsigned int i,
                             const arma::vec &p,
                             unsigned int comp) const
{
    ASSERT_DBG(i < dim_).error("Index of basis function is out of range.");
    ASSERT_DBG(comp < n_components()).error("Index of component is out of range.");

    // component index in the base FE
    int l_comp = comp-dof_indices_[i].component_offset;
    if (l_comp >= 0 && l_comp < static_cast<int>(fs_[dof_indices_[i].fe_index]->n_components()))
        return fs_[dof_indices_[i].fe_index]->basis_value(dof_indices_[i].basis_index, p, l_comp);
    else
        return 0;
}


const arma::vec FESystemFunctionSpace::basis_grad(const unsigned int i, 
                                                  const arma::vec &p, 
                                                  const unsigned int comp) const
{
    ASSERT_DBG(i < dim_).error("Index of basis function is out of range.");
    ASSERT_DBG(comp < n_components()).error("Index of component is out of range.");

    // component index in the base FE
    int l_comp = comp-dof_indices_[i].component_offset;
    if (l_comp >= 0 && l_comp < static_cast<int>(fs_[dof_indices_[i].fe_index]->n_components()))
        return fs_[dof_indices_[i].fe_index]->basis_grad(dof_indices_[i].basis_index, p, l_comp);
    else
        return arma::zeros(space_dim());
}




template<unsigned int dim>
FESystem<dim>::FESystem(std::shared_ptr<FiniteElement<dim> > fe, FEType t)
{
  OLD_ASSERT(fe->n_components() == 1, "FEVectorContravariant and FEVectorPiola can only be created from scalar FE.");
  OLD_ASSERT(t == FEType::FEVectorContravariant || t == FEType::FEVectorPiola,
             "This constructor can be used only for FEVectorContravariant or FEVectorPiola.");
  
  FiniteElement<dim>::init(false, t);
  fe_ = std::vector<std::shared_ptr<FiniteElement<dim> > >(dim, fe);
  initialize();
}


template<unsigned int dim>
FESystem<dim>::FESystem(const std::shared_ptr<FiniteElement<dim> > &fe, FEType t, unsigned int n)
{
    OLD_ASSERT(t == FEType::FEVector || t == FEType::FETensor || t == FEType::FEMixedSystem,
               "This constructor can be used only for FEVector, FETensor or FEMixedSystem.");
    OLD_ASSERT(fe->n_components() == 1 || t == FEType::FEMixedSystem,
               "FEVector and FETensor can only be created from scalar FE.");
    
    FiniteElement<dim>::init(false, t);
    fe_ = std::vector<std::shared_ptr<FiniteElement<dim> > >(n, fe);
    initialize();
}


template<unsigned int dim>
FESystem<dim>::FESystem(std::vector<std::shared_ptr<FiniteElement<dim> > > fe)
{
  FiniteElement<dim>::init(false, FEMixedSystem);
  for (std::shared_ptr<FiniteElement<dim> > fe_object : fe)
    fe_.push_back(fe_object);
  initialize();
}



template<unsigned int dim>
void FESystem<dim>::initialize()
{
  unsigned int fe_index = 0;
  unsigned int comp_offset = 0;
  vector<shared_ptr<FunctionSpace> > fs_vector;
  // for each base FE add components, support points, and other 
  // information to the system
  for (auto fe : fe_)
  {
    switch (fe->type_)
    {
      case FEType::FEScalar:
        scalar_components_.push_back(comp_offset);
        break;
      case FEType::FEVector:
      case FEType::FEVectorContravariant:
      case FEType::FEVectorPiola:
        vector_components_.push_back(comp_offset);
        break;
      case FEType::FETensor:
        tensor_components_.push_back(comp_offset);
        break;
      default:
        OLD_ASSERT(false, "Not implemented.");
        break;
    }

    fe_index++;
    comp_offset += fe->n_components();
    fs_vector.push_back(shared_ptr<FunctionSpace>(fe->function_space_));
  }
  
  this->function_space_ = make_shared<FESystemFunctionSpace>(fs_vector);
  
  double dof_index = 0;
  comp_offset = 0;
  for (auto fe : fe_)
  {
      for (unsigned int i=0; i<fe->n_dofs(); i++)
      {
          arma::vec coefs(this->function_space_->n_components());
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
      for (uint i=0; i<fe->n_dofs(); ++i)
      {
        std::vector<bool> nonzeros(this->function_space_->n_components(), false);
        for (unsigned int c=0; c<fe->n_components(); c++)\
          nonzeros[comp_offset+c] = fe->get_nonzero_components(i)[c];
        this->nonzero_components_.push_back(nonzeros);
        dof_index++;
      }
      comp_offset += fe->n_components();
    }
//   }

  compute_node_matrix();
}




template<unsigned int dim> inline
UpdateFlags FESystem<dim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    for (auto fe : fe_)
      f |= fe->update_each(flags);

    return f;
}

template<unsigned int dim>
vector< arma::vec::fixed<dim+1> > FESystem<dim>::dof_points() const {
    std::vector<arma::vec::fixed<dim+1>> points(20);
    points.resize(0);
    for (unsigned int i = 0; i < fe_.size(); i++)
    {
        auto fe_points = fe_[i]->dof_points();
        points.insert(points.end(), fe_points.begin(), fe_points.end());
    }
    return points;
}


template<unsigned int dim>
void FESystem<dim>::compute_node_matrix()
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




template<unsigned int dim>
std::vector<unsigned int> FESystem<dim>::fe_dofs(unsigned int fe_index)
{
    auto fs = dynamic_cast<FESystemFunctionSpace *>(this->function_space_.get());
    std::vector<unsigned int> fe_dof_indices;
    
    for (unsigned int i=0; i<fs->dof_indices().size(); i++)
    {
        if (fs->dof_indices()[i].fe_index == fe_index)
            fe_dof_indices.push_back(i);
    }
    
    return fe_dof_indices;
}








template class FESystem<0>;
template class FESystem<1>;
template class FESystem<2>;
template class FESystem<3>;



