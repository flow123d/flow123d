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
 * @file    fe_values.cc
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel
 */

#include "fem/mapping_p1.hh"
#include "quadrature/quadrature.hh"
#include "fem/element_values.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"
#include "fem/fe_system.hh"



using namespace arma;
using namespace std;







template<unsigned int spacedim>
FEValues<spacedim>::FEInternalData::FEInternalData(unsigned int np, unsigned int nd)
    : n_points(np),
      n_dofs(nd)
{
    ref_shape_values.resize(np, vector<arma::vec>(nd));
    ref_shape_grads.resize(np, vector<arma::mat>(nd));
}


template<unsigned int spacedim>
FEValues<spacedim>::FEInternalData::FEInternalData(const FEInternalData &fe_system_data,
                               const std::vector<unsigned int> &dof_indices,
                               unsigned int first_component_idx,
                               unsigned int ncomps)
    : FEInternalData(fe_system_data.n_points, dof_indices.size())
{
    for (unsigned int ip=0; ip<n_points; ip++)
        for (unsigned int id=0; id<dof_indices.size(); id++)
        {
            ref_shape_values[ip][id] = fe_system_data.ref_shape_values[ip][dof_indices[id]].subvec(first_component_idx, first_component_idx+ncomps-1);
            ref_shape_grads[ip][id] = fe_system_data.ref_shape_grads[ip][dof_indices[id]].cols(first_component_idx, first_component_idx+ncomps-1);
        }
}



template<unsigned int spacedim>
template<unsigned int DIM>
void FEValues<spacedim>::ViewsCache::initialize(const FEValues<spacedim> &fv, const FiniteElement<DIM> &fe)
{
  scalars.clear();
  vectors.clear();
  tensors.clear();
  switch (fe.type_) {
    case FEType::FEScalar:
      scalars.push_back(FEValuesViews::Scalar<spacedim>(fv, 0));
      break;
    case FEType::FEVector:
    case FEType::FEVectorContravariant:
    case FEType::FEVectorPiola:
      vectors.push_back(FEValuesViews::Vector<spacedim>(fv, 0));
      break;
    case FEType::FETensor:
      tensors.push_back(FEValuesViews::Tensor<spacedim>(fv, 0));
      break;
    case FEType::FEMixedSystem:
      const FESystem<DIM> *fe_sys = dynamic_cast<const FESystem<DIM>*>(&fe);
      ASSERT_DBG(fe_sys != nullptr).error("Mixed system must be represented by FESystem.");
      
      // Loop through sub-elements and add views according to their types.
      // Note that the component index is calculated using fe->n_space_components(),
      // not fe->n_components()!
      unsigned int comp_offset = 0;
      for (auto fe : fe_sys->fe())
      {
          switch (fe->type_)
          {
          case FEType::FEScalar:
              scalars.push_back(FEValuesViews::Scalar<spacedim>(fv,comp_offset));
              break;
          case FEType::FEVector:
          case FEType::FEVectorContravariant:
          case FEType::FEVectorPiola:
              vectors.push_back(FEValuesViews::Vector<spacedim>(fv,comp_offset));
              break;
          case FEType::FETensor:
              tensors.push_back(FEValuesViews::Tensor<spacedim>(fv,comp_offset));
              break;
          default:
              ASSERT_DBG(false).error("Not implemented.");
              break;
          }

          comp_offset += fe->n_space_components(spacedim);
      }
      break;
  }
}



template<unsigned int spacedim>
FEValues<spacedim>::FEValues()
: dim_(-1), n_points_(0), n_dofs_(0), elm_values(nullptr), fe_data(nullptr)
{
}



template<unsigned int spacedim>
FEValues<spacedim>::~FEValues() {
    if (elm_values != nullptr) delete elm_values;
    if (fe_data) delete fe_data;
    for (unsigned int sid=0; sid<side_fe_data.size(); sid++)
        for (unsigned int pid=0; pid<side_fe_data[sid].size(); pid++)
            delete side_fe_data[sid][pid];
}



template<unsigned int spacedim>
template<unsigned int DIM>
void FEValues<spacedim>::allocate(
        unsigned int n_points,
        FiniteElement<DIM> & _fe,
        UpdateFlags _flags)
{
    // For FEVector and FETensor check number of components.
    // This cannot be done in FiniteElement since it does not know spacedim.
    if (_fe.type_ == FEVector) {
        ASSERT_DBG(_fe.n_components() == spacedim).error("FEVector must have spacedim components.");
    } else if (_fe.type_ == FETensor) {
        ASSERT_DBG(_fe.n_components() == spacedim*spacedim).error("FETensor must have spacedim*spacedim components.");
    }
    
    dim_ = DIM;
    n_points_ = n_points;
    n_dofs_ = _fe.n_dofs();
    n_components_ = _fe.n_space_components(spacedim);
    fe_type_ = _fe.type_;
    FESystem<DIM> *fe_sys = dynamic_cast<FESystem<DIM>*>(&_fe);
    if (fe_sys != nullptr)
    {
        for (unsigned int f=0; f<fe_sys->fe().size(); f++)
        {
            fe_sys_dofs_.push_back(fe_sys->fe_dofs(f));
            fe_sys_n_components_.push_back(fe_sys->fe()[f]->n_components());
            fe_sys_n_space_components_.push_back(fe_sys->fe()[f]->n_space_components(spacedim));
        }
    }
    
    // add flags required by the finite element or mapping
    update_flags = _flags | _fe.update_each(_flags);
    update_flags |= MappingP1<DIM,spacedim>::update_each(update_flags);
    if (update_flags & update_values)
        shape_values.resize(n_points_, vector<double>(n_dofs_*n_components_));

    if (update_flags & update_gradients)
        shape_gradients.resize(n_points_, vector<arma::vec::fixed<spacedim> >(n_dofs_*n_components_));
    
    views_cache_.initialize(*this, _fe);
}



template<unsigned int spacedim>
template<unsigned int DIM>
typename FEValues<spacedim>::FEInternalData *FEValues<spacedim>::init_fe_data(const FiniteElement<DIM> &fe, const Quadrature &q)
{
    ASSERT_DBG( DIM == dim_ );
    ASSERT_DBG( q.dim() == DIM );
    FEInternalData *data = new FEInternalData(q.size(), n_dofs_);

    arma::mat shape_values(n_dofs_, fe.n_components());
    for (unsigned int i=0; i<q.size(); i++)
    {
        for (unsigned int j=0; j<n_dofs_; j++)
        {
            for (unsigned int c=0; c<fe.n_components(); c++)
                shape_values(j,c) = fe.shape_value(j, q.point<DIM>(i), c);
            
            data->ref_shape_values[i][j] = trans(shape_values.row(j));
        }
    }

    arma::mat grad(DIM, fe.n_components());
    for (unsigned int i=0; i<q.size(); i++)
    {
        for (unsigned int j=0; j<n_dofs_; j++)
        {
            grad.zeros();
            for (unsigned int c=0; c<fe.n_components(); c++)
                grad.col(c) += fe.shape_grad(j, q.point<DIM>(i), c);
            
            data->ref_shape_grads[i][j] = grad;
        }
    }
    
    return data;
}


template<unsigned int spacedim>
double FEValues<spacedim>::shape_value(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, n_points_);
  return shape_values[point_no][function_no];
}


template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValues<spacedim>::shape_grad(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, n_points_);
  return shape_gradients[point_no][function_no];
}


template<unsigned int spacedim>
double FEValues<spacedim>::shape_value_component(const unsigned int function_no, 
                                    const unsigned int point_no, 
                                    const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, n_points_);
  ASSERT_LT_DBG(comp, n_components_);
  return shape_values[point_no][function_no*n_components_+comp];
}


template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValues<spacedim>::shape_grad_component(const unsigned int function_no,
                                                        const unsigned int point_no,
                                                        const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, n_points_);
  ASSERT_LT_DBG(comp, n_components_);
  return shape_gradients[point_no][function_no*n_components_+comp];
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_scalar_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FEScalar);
    
    // shape values
    if (update_flags & update_values)
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                shape_values[i][j] = fe_data.ref_shape_values[i][j][0];

    // shape gradients
    if (update_flags & update_gradients)
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                shape_gradients[i][j] = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_vec_data(const ElementValues<spacedim> &elm_values,
                                           const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FEVector);
    
    // shape values
    if (update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_vec_contravariant_data(const ElementValues<spacedim> &elm_values,
                                                         const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FEVectorContravariant);
    
    // shape values
    if (update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = elm_values.jacobian(i) * fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j] * trans(elm_values.jacobian(i));
                for (unsigned int c=0; c<spacedim; c++)
                    shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_vec_piola_data(const ElementValues<spacedim> &elm_values,
                                                 const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FEVectorPiola);
    
    // shape values
    if (update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = elm_values.jacobian(i)*fe_data.ref_shape_values[i][j]/elm_values.determinant(i);
                for (unsigned int c=0; c<spacedim; c++)
                    shape_values[i][j*spacedim+c] = fv_vec(c);
            }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j] * trans(elm_values.jacobian(i))
                        / elm_values.determinant(i);
                for (unsigned int c=0; c<spacedim; c++)
                    shape_gradients[i][j*spacedim+c] = grads.col(c);
            }   
    }
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_tensor_data(const ElementValues<spacedim> &elm_values,
                                              const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FETensor);
    
    // shape values
    if (update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    shape_values[i][j*spacedim*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    shape_gradients[i][j*spacedim*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_system_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data)
{
    ASSERT_DBG(fe_type_ == FEMixedSystem);
    
    // for mixed system we first fill data in sub-elements
    unsigned int comp_offset = 0;
    for (unsigned int f=0; f<fe_sys_dofs_.size(); f++)
    {
        // fill fe_values for base FE
        FEInternalData vec_fe_data(fe_data, fe_sys_dofs_[f], comp_offset, fe_sys_n_components_[f]);
        fe_values_vec[f].fill_data(elm_values, vec_fe_data);
        
        comp_offset += fe_sys_n_components_[f];
    }
    
    // shape values
    if (update_flags & update_values)
    {
        arma::vec fv_vec;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys_dofs_.size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys_dofs_[f].size(); n++)
                    for (unsigned int c=0; c<fe_sys_n_space_components_[f]; c++)
                        shape_values[i][shape_offset+n_components_*n+comp_offset+c] = fe_values_vec[f].shape_values[i][n*fe_sys_n_space_components_[f]+c];
            
            comp_offset += fe_sys_n_space_components_[f];
            shape_offset += fe_sys_dofs_[f].size()*n_components_;
        }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        arma::mat grads;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys_dofs_.size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys_dofs_[f].size(); n++)
                    for (unsigned int c=0; c<fe_sys_n_space_components_[f]; c++)
                        shape_gradients[i][shape_offset+n_components_*n+comp_offset+c] = fe_values_vec[f].shape_gradients[i][n*fe_sys_n_space_components_[f]+c];
            
            comp_offset += fe_sys_n_space_components_[f];
            shape_offset += fe_sys_dofs_[f].size()*n_components_;
        }
    }
    
}


template<unsigned int spacedim>
void FEValues<spacedim>::fill_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data)
{
    switch (fe_type_) {
        case FEScalar:
            fill_scalar_data(elm_values, fe_data);
            break;
        case FEVector:
            fill_vec_data(elm_values, fe_data);
            break;
        case FEVectorContravariant:
            fill_vec_contravariant_data(elm_values, fe_data);
            break;
        case FEVectorPiola:
            fill_vec_piola_data(elm_values, fe_data);
            break;
        case FETensor:
            fill_tensor_data(elm_values, fe_data);
            break;
        case FEMixedSystem:
            fill_system_data(elm_values, fe_data);
            break;
        default:
            ASSERT(false).error("Not implemented.");
    }
}








template<unsigned int spacedim>
template<unsigned int DIM>
void FEValues<spacedim>::initialize(
         Quadrature &q,
         FiniteElement<DIM> &_fe,
         UpdateFlags _flags)
{
    if (DIM == 0) return; // avoid unnecessary allocation of dummy 0 dimensional objects

    allocate( q.size(), _fe, _flags);
    elm_values = new ElementValues<spacedim>(q, update_flags, DIM);

    // In case of mixed system allocate data for sub-elements.
    if (fe_type_ == FEMixedSystem)
    {
        FESystem<DIM> *fe = dynamic_cast<FESystem<DIM>*>(&_fe);
        ASSERT_DBG(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        fe_values_vec.resize(fe->fe().size());
        for (unsigned int f=0; f<fe->fe().size(); f++)
            fe_values_vec[f].initialize(q, *fe->fe()[f], update_flags);
    }

    // precompute finite element data
    if ( q.dim() == DIM )
    {
        fe_data = init_fe_data(_fe, q);
    }
    else if ( q.dim() + 1 == DIM )
    {
        side_fe_data.resize(RefElement<DIM>::n_sides);
        for (unsigned int sid = 0; sid < RefElement<DIM>::n_sides; sid++)
        {
            side_fe_data[sid].resize(RefElement<DIM>::n_side_permutations);

            // For each side transform the side quadrature points to the cell quadrature points
            // and then precompute side_fe_data.
            for (unsigned int pid = 0; pid < RefElement<DIM>::n_side_permutations; pid++)
                side_fe_data[sid][pid] = init_fe_data(_fe, q.make_from_side<DIM>(sid,pid));
        }
    }
    else
        ASSERT_DBG(false)(q.dim())(DIM).error("Dimension mismatch in FEValues::initialize().");
}





template<unsigned int spacedim>
void FEValues<spacedim>::reinit(const ElementAccessor<spacedim> &cell)
{
	ASSERT_EQ_DBG( dim_, cell.dim() );
    
    if (!elm_values->cell().is_valid() ||
        elm_values->cell() != cell)
    {
        elm_values->reinit(cell);
    }
    
    fill_data(*elm_values, *fe_data);
}


template<unsigned int spacedim>
void FEValues<spacedim>::reinit(const Side &cell_side)
{
    ASSERT_EQ_DBG( dim_, cell_side.dim() );
    
    if (!elm_values->side().is_valid() || 
        elm_values->side() != cell_side)
    {
        elm_values->reinit(cell_side);
    }

    const LongIdx sid = cell_side.side_idx();
    const unsigned int pid = elm_values->side().element()->permutation_idx(sid);
    
    // calculation of finite element data
    fill_data(*elm_values, *side_fe_data[sid][pid]);
}



std::vector<FEValues<3>> mixed_fe_values(
        QGauss::array &quadrature,
        MixedPtr<FiniteElement> fe,
        UpdateFlags flags)
{
    std::vector<FEValues<3>> fv(4);
    fv[0].initialize(quadrature[0], *fe.get<0>(), flags);
    fv[1].initialize(quadrature[1], *fe.get<1>(), flags);
    fv[2].initialize(quadrature[2], *fe.get<2>(), flags);
    fv[3].initialize(quadrature[3], *fe.get<3>(), flags);
    return fv;
}




















template class FEValues<3>;
