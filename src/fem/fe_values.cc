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
#include "fem/fe_values_map.hh"



using namespace arma;
using namespace std;







template<class FV, unsigned int spacedim>
FEValuesBase<FV, spacedim>::FEValuesBase()
: dim_(-1), n_points_(0), n_dofs_(0)
{
}


FEInternalData::FEInternalData(unsigned int np, unsigned int nd)
    : n_points(np),
      n_dofs(nd)
{
    ref_shape_values.resize(np, vector<arma::vec>(nd));
    ref_shape_grads.resize(np, vector<arma::mat>(nd));
}


FEInternalData::FEInternalData(const FEInternalData &fe_system_data,
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



template<class FV, unsigned int spacedim>
template<unsigned int DIM>
void FEValuesBase<FV, spacedim>::ViewsCache::initialize(const FV &fv, const FiniteElement<DIM> &fe)
{
  scalars.clear();
  vectors.clear();
  tensors.clear();
  switch (fe.type_) {
    case FEType::FEScalar:
      scalars.push_back(FEValuesViews::Scalar<FV, spacedim>(fv, 0));
      break;
    case FEType::FEVector:
    case FEType::FEVectorContravariant:
    case FEType::FEVectorPiola:
      vectors.push_back(FEValuesViews::Vector<FV, spacedim>(fv, 0));
      break;
    case FEType::FETensor:
      tensors.push_back(FEValuesViews::Tensor<FV, spacedim>(fv, 0));
      break;
    case FEType::FEMixedSystem:
      const FESystem<DIM> *fe_sys = dynamic_cast<const FESystem<DIM>*>(&fe);
      ASSERT(fe_sys != nullptr).error("Mixed system must be represented by FESystem.");
      
      // Loop through sub-elements and add views according to their types.
      // Note that the component index is calculated using fe->n_space_components(),
      // not fe->n_components()!
      unsigned int comp_offset = 0;
      for (auto fe : fe_sys->fe())
      {
          switch (fe->type_)
          {
          case FEType::FEScalar:
              scalars.push_back(FEValuesViews::Scalar<FV, spacedim>(fv,comp_offset));
              break;
          case FEType::FEVector:
          case FEType::FEVectorContravariant:
          case FEType::FEVectorPiola:
              vectors.push_back(FEValuesViews::Vector<FV, spacedim>(fv,comp_offset));
              break;
          case FEType::FETensor:
              tensors.push_back(FEValuesViews::Tensor<FV, spacedim>(fv,comp_offset));
              break;
          default:
              ASSERT(false).error("Not implemented.");
              break;
          }

          comp_offset += fe->n_space_components(spacedim);
      }
      break;
  }
}



template<class FV, unsigned int spacedim>
template<unsigned int DIM>
void FEValuesBase<FV, spacedim>::initialize(
         Quadrature &q,
         FiniteElement<DIM> &_fe,
         UpdateFlags _flags)
{
    if (DIM == 0) //return; // avoid unnecessary allocation of dummy 0 dimensional objects
    	ASSERT(q.size() == 1);

    this->allocate( q, _fe, _flags);
    this->initialize_in(q, DIM);

    // In case of mixed system allocate data for sub-elements.
    if (fe_type_ == FEMixedSystem)
    {
        FESystem<DIM> *fe = dynamic_cast<FESystem<DIM>*>(&_fe);
        ASSERT(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        fe_values_vec.resize(fe->fe().size());
        init_fe_val_vec();
        for (unsigned int f=0; f<fe->fe().size(); f++)
            fe_values_vec[f].initialize(q, *fe->fe()[f], update_flags);
    }

    // precompute finite element data
    if ( q.dim() == DIM )
    {
        fe_data_ = init_fe_data(_fe, q);
    }
    else if ( q.dim() + 1 == DIM )
    {
        side_fe_data_.resize(RefElement<DIM>::n_sides);
        for (unsigned int sid = 0; sid < RefElement<DIM>::n_sides; sid++)
        {
            side_fe_data_[sid] = init_fe_data(_fe, q.make_from_side<DIM>(sid));
        }
    }
    else
        ASSERT(false)(q.dim())(DIM).error("Dimension mismatch in FEValues::initialize().");
}



template<class FV, unsigned int spacedim>
template<unsigned int DIM>
void FEValuesBase<FV, spacedim>::allocate(
        Quadrature &_q,
        FiniteElement<DIM> & _fe,
        UpdateFlags _flags)
{
    // For FEVector and FETensor check number of components.
    // This cannot be done in FiniteElement since it does not know spacedim.
    if (_fe.type_ == FEVector) {
        ASSERT(_fe.n_components() == spacedim).error("FEVector must have spacedim components.");
    } else if (_fe.type_ == FETensor) {
        ASSERT(_fe.n_components() == spacedim*spacedim).error("FETensor must have spacedim*spacedim components.");
    }

    fe_sys_dofs_.clear();
    fe_sys_n_components_.clear();
    fe_sys_n_space_components_.clear();
    
    dim_ = DIM;
    n_points_ = _q.size();
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
    this->allocate_in(_q.dim());

    views_cache_.initialize(*this->fv_, _fe);
}



template<class FV, unsigned int spacedim>
template<unsigned int DIM>
std::shared_ptr<FEInternalData> FEValuesBase<FV, spacedim>::init_fe_data(const FiniteElement<DIM> &fe, const Quadrature &q)
{
    ASSERT( DIM == dim_ );
    ASSERT( q.dim() == DIM );
    std::shared_ptr<FEInternalData> data = std::make_shared<FEInternalData>(q.size(), n_dofs_);

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


template<class FV, unsigned int spacedim>
void FEValuesBase<FV, spacedim>::fill_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data)
{
    switch (this->fe_type_) {
        case FEScalar:
            this->fill_data_specialized<MapScalar<FV, spacedim>>(elm_values, fe_data);
            break;
        case FEVector:
            this->fill_data_specialized<MapVector<FV, spacedim>>(elm_values, fe_data);
            break;
        case FEVectorContravariant:
            this->fill_data_specialized<MapContravariant<FV, spacedim>>(elm_values, fe_data);
            break;
        case FEVectorPiola:
            this->fill_data_specialized<MapPiola<FV, spacedim>>(elm_values, fe_data);
            break;
        case FETensor:
            this->fill_data_specialized<MapTensor<FV, spacedim>>(elm_values, fe_data);
            break;
        case FEMixedSystem:
            this->fill_data_specialized<MapSystem<FV, spacedim>>(elm_values, fe_data);
            break;
        default:
            ASSERT_PERMANENT(false).error("Not implemented.");
    }
}



template<class FV, unsigned int spacedim>
template<class MapType>
inline void FEValuesBase<FV, spacedim>::fill_data_specialized(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data) {
	MapType map_type;
	map_type.fill_values_vec(*this->fv_, elm_values, fe_data);
    if (this->update_flags & update_values)
    	map_type.update_values(*this->fv_, elm_values, fe_data);
    if (this->update_flags & update_gradients)
    	map_type.update_gradients(*this->fv_, elm_values, fe_data);
}




/*template<unsigned int spacedim>
double FEValues<spacedim>::shape_value_component(const unsigned int function_no, 
                                    const unsigned int point_no, 
                                    const unsigned int comp) const
{
  ASSERT_LT(function_no, n_dofs_);
  ASSERT_LT(point_no, n_points_);
  ASSERT_LT(comp, n_components_);
  return shape_values[point_no][function_no*n_components_+comp];
}*/


template<unsigned int spacedim>
FEValues<spacedim>::FEValues()
: FEValuesBase<FEValues<spacedim>, spacedim>() {}


template<unsigned int spacedim>
FEValues<spacedim>::~FEValues() {
}


template<unsigned int spacedim>
void FEValues<spacedim>::initialize_in (
         Quadrature &q,
		 unsigned int dim)
{
    elm_values_ = std::make_shared<ElementValues<spacedim> >(q, this->update_flags, dim);
}



template<unsigned int spacedim>
void FEValues<spacedim>::allocate_in(FMT_UNUSED unsigned int q_dim)
{
    if (this->update_flags & update_values)
        shape_values_.resize(this->n_points_, vector<double>(this->n_dofs_*this->n_components_));

    if (this->update_flags & update_gradients)
        shape_gradients_.resize(this->n_points_, vector<arma::vec::fixed<spacedim> >(this->n_dofs_*this->n_components_));

    this->fv_ = this;
}



template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValues<spacedim>::shape_grad_component(const unsigned int function_no,
                                                        const unsigned int point_no,
                                                        const unsigned int comp) const
{
  ASSERT_LT(function_no, this->n_dofs_);
  ASSERT_LT(point_no, this->n_points_);
  ASSERT_LT(comp, this->n_components_);
  return shape_gradients_[point_no][function_no*this->n_components_+comp];
}


/*template<unsigned int spacedim>
void FEValues<spacedim>::fill_scalar_data(const ElementValues<spacedim> &elm_values, const FEInternalData &fe_data)
{
    ASSERT(fe_type_ == FEScalar);
    
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
    ASSERT(fe_type_ == FEVector);
    
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
    ASSERT(fe_type_ == FEVectorContravariant);
    
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
    ASSERT(fe_type_ == FEVectorPiola);
    
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
    ASSERT(fe_type_ == FETensor);
    
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
    ASSERT(fe_type_ == FEMixedSystem);
    
    // for mixed system we first fill data in sub-elements
    unsigned int comp_offset = 0;
    for (unsigned int f=0; f<fe_sys_dofs_.size(); f++)
    {
        // fill fe_values for base FE
        FEInternalData vec_fe_data(fe_data, fe_sys_dofs_[f], comp_offset, fe_sys_n_components_[f]);
        fe_values_vec[f].fill_data(elm_values, vec_fe_data); // fe_values.fe_values_vec
        
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
    
}*/


template<unsigned int spacedim>
void FEValues<spacedim>::reinit(const ElementAccessor<spacedim> &cell)
{
	ASSERT_EQ( this->dim_, cell.dim() );
    
    if (!elm_values_->cell().is_valid() ||
        elm_values_->cell() != cell)
    {
        elm_values_->reinit(cell);
    }
    
    this->fill_data(*elm_values_, *this->fe_data_);
}


template<unsigned int spacedim>
void FEValues<spacedim>::reinit(const Side &cell_side)
{
    ASSERT_EQ( this->dim_, cell_side.dim()+1 );
    
    if (!elm_values_->side().is_valid() ||
        elm_values_->side() != cell_side)
    {
        elm_values_->reinit(cell_side);
    }

    const LongIdx sid = cell_side.side_idx();
    
    // calculation of finite element data
    this->fill_data(*elm_values_, *(this->side_fe_data_[sid]) );
}



//template<unsigned int spacedim>
//PatchFEValues_TEMP<spacedim>::PatchFEValues_TEMP(unsigned int max_size)
//: FEValuesBase<PatchFEValues_TEMP<spacedim>, spacedim>(),
//  patch_data_idx_(-1), used_size_(0), max_n_elem_(max_size) {}
//
//
//template<unsigned int spacedim>
//void PatchFEValues_TEMP<spacedim>::reinit(PatchElementsList patch_elements) {
//    element_patch_map_.clear();
//    if (object_type_ == ElementFE)
//        used_size_ = patch_elements.size();
//    else
//        used_size_ = patch_elements.size() * (this->dim_+1);
//    ASSERT_LE(used_size_, max_size());
//
//    unsigned int i=0;
//    for (auto it=patch_elements.begin(); it!=patch_elements.end(); ++it, ++i) {
//        if (object_type_ == ElementFE) {
//            patch_data_idx_ = i;
//            element_patch_map_[it->second] = i;
//            element_data_[i].elm_values_->reinit(it->first);
//            this->fill_data(*element_data_[i].elm_values_, *this->fe_data_);
//        } else {
//            element_patch_map_[it->second] = i * (this->dim_+1);
//            for (unsigned int sid=0; sid<this->dim_+1; ++sid) {
//                patch_data_idx_ = i * (this->dim_+1) + sid;
//                element_data_[patch_data_idx_].elm_values_->reinit( *it->first.side(sid) );
//                this->fill_data(*element_data_[patch_data_idx_].elm_values_, *this->side_fe_data_[sid]);
//
//            }
//        }
//    }
//}
//
//
//template<unsigned int spacedim>
//void PatchFEValues_TEMP<spacedim>::allocate_in(unsigned int q_dim)
//{
//    ASSERT_PERMANENT_GT(this->max_n_elem_, 0);
//
//    if ( q_dim == this->dim_ ) {
//        element_data_.resize( this->max_n_elem_ );
//        object_type_ = ElementFE;
//    } else if ( q_dim+1 == this->dim_ ) {
//        element_data_.resize( this->max_n_elem_ * (this->dim_+1) );
//        object_type_ = SideFE;
//    } else
//        ASSERT(false)(q_dim)(this->dim_).error("Invalid dimension of quadrature!");
//
//    for (uint i=0; i<max_size(); ++i) {
//        if (this->update_flags & update_values)
//            element_data_[i].shape_values_.resize(this->n_points_, vector<double>(this->n_dofs_*this->n_components_));
//
//        if (this->update_flags & update_gradients)
//            element_data_[i].shape_gradients_.resize(this->n_points_, vector<arma::vec::fixed<spacedim> >(this->n_dofs_*this->n_components_));
//    }
//
//    this->fv_ = this;
//}
//
//
//template<unsigned int spacedim>
//void PatchFEValues_TEMP<spacedim>::initialize_in(
//        Quadrature &q,
//        unsigned int dim)
//{
//    for (uint i=0; i<max_size(); ++i)
//        element_data_[i].elm_values_ = std::make_shared<ElementValues<spacedim> >(q, this->update_flags, dim);
//}
//
//
//template<unsigned int spacedim>
//void PatchFEValues_TEMP<spacedim>::init_fe_val_vec()
//{
//    for (unsigned int i=0; i<this->fe_values_vec.size(); ++i)
//        this->fe_values_vec[i].resize( this->max_size() );
//}
//
//
//template<unsigned int spacedim>
//arma::vec::fixed<spacedim> PatchFEValues_TEMP<spacedim>::shape_grad_component(const unsigned int function_no,
//                                                        const unsigned int point_no,
//                                                        const unsigned int comp) const
//{
//    ASSERT_LT(function_no, this->n_dofs_);
//    ASSERT_LT(point_no, this->n_points_);
//    ASSERT_LT(comp, this->n_components_);
//    return element_data_[patch_data_idx_].shape_gradients_[point_no][function_no*this->n_components_+comp];
//}



std::vector<FEValues<3>> mixed_fe_values(
        QGauss::array &quadrature,
        MixedPtr<FiniteElement> fe,
        UpdateFlags flags)
{
    std::vector<FEValues<3>> fv(4);
    fv[0].initialize(quadrature[0], *fe[0_d], flags);
    fv[1].initialize(quadrature[1], *fe[1_d], flags);
    fv[2].initialize(quadrature[2], *fe[2_d], flags);
    fv[3].initialize(quadrature[3], *fe[3_d], flags);
    return fv;
}



// explicit instantiation
template void FEValuesBase<FEValues<3>, 3>::initialize<0>(Quadrature&, FiniteElement<0>&, UpdateFlags);
template void FEValuesBase<FEValues<3>, 3>::initialize<1>(Quadrature&, FiniteElement<1>&, UpdateFlags);
template void FEValuesBase<FEValues<3>, 3>::initialize<2>(Quadrature&, FiniteElement<2>&, UpdateFlags);
template void FEValuesBase<FEValues<3>, 3>::initialize<3>(Quadrature&, FiniteElement<3>&, UpdateFlags);

template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapScalar<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);
template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapPiola<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);
template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapContravariant<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);
template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapVector<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);
template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapTensor<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);
template void FEValuesBase<FEValues<3>, 3>::fill_data_specialized<MapSystem<FEValues<3>, 3>>(const ElementValues<3> &, const FEInternalData &);






template class FEValues<3>;
