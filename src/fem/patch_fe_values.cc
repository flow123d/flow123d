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
 * @file    patch_fe_values.cc
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel, David Flanderka
 */

#include "fem/patch_fe_values.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"



//template<unsigned int spacedim>
//PatchFEValues<spacedim>::FEInternalData::FEInternalData(unsigned int np, unsigned int nd)
//    : n_points(np),
//      n_dofs(nd)
//{
//    ref_shape_values.resize(np, vector<arma::vec>(nd));
//    ref_shape_grads.resize(np, vector<arma::mat>(nd));
//}
//
//
//template<unsigned int spacedim>
//PatchFEValues<spacedim>::FEInternalData::FEInternalData(const PatchFEValues<spacedim>::FEInternalData &fe_system_data,
//                               const std::vector<unsigned int> &dof_indices,
//                               unsigned int first_component_idx,
//                               unsigned int ncomps)
//    : FEInternalData(fe_system_data.n_points, dof_indices.size())
//{
//    for (unsigned int ip=0; ip<n_points; ip++)
//        for (unsigned int id=0; id<dof_indices.size(); id++)
//        {
//            ref_shape_values[ip][id] = fe_system_data.ref_shape_values[ip][dof_indices[id]].subvec(first_component_idx, first_component_idx+ncomps-1);
//            ref_shape_grads[ip][id] = fe_system_data.ref_shape_grads[ip][dof_indices[id]].cols(first_component_idx, first_component_idx+ncomps-1);
//        }
//}


template<unsigned int spacedim>
template<unsigned int DIM>
void PatchFEValues<spacedim>::DimPatchFEValues::initialize(
         Quadrature &q,
         FiniteElement<DIM> &_fe,
         UpdateFlags _flags)
{
    if (DIM == 0) //return; // avoid unnecessary allocation of dummy 0 dimensional objects
    	ASSERT(q.size() == 1);

    this->allocate( q, _fe, _flags);
    for (uint i=0; i<max_size(); ++i)
        element_data_[i].elm_values_ = std::make_shared<ElementValues<spacedim> >(q, this->update_flags, DIM);

    // In case of mixed system allocate data for sub-elements.
//    if (fe_type_ == FEMixedSystem)
//    {
//        FESystem<DIM> *fe = dynamic_cast<FESystem<DIM>*>(&_fe);
//        ASSERT(fe != nullptr).error("Mixed system must be represented by FESystem.");
//
//        fe_values_vec.resize(fe->fe().size());
//        init_fe_val_vec();
//        for (unsigned int f=0; f<fe->fe().size(); f++)
//            fe_values_vec[f].initialize(q, *fe->fe()[f], update_flags);
//    }

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


template<unsigned int spacedim>
template<unsigned int DIM>
void PatchFEValues<spacedim>::DimPatchFEValues::allocate(
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
    ASSERT_PERMANENT_GT(this->max_n_elem_, 0);

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

    if ( _q.dim() == this->dim_ ) {
        element_data_.resize( this->max_n_elem_ );
        object_type_ = ElementFE;
    } else if ( _q.dim()+1 == this->dim_ ) {
        element_data_.resize( this->max_n_elem_ * (this->dim_+1) );
        object_type_ = SideFE;
    } else
        ASSERT(false)(_q.dim())(this->dim_).error("Invalid dimension of quadrature!");

    for (uint i=0; i<max_size(); ++i) {
        if (this->update_flags & update_values)
            element_data_[i].shape_values_.resize(this->n_points_, vector<double>(this->n_dofs_*this->n_components_));

        if (this->update_flags & update_gradients)
            element_data_[i].shape_gradients_.resize(this->n_points_, vector<arma::vec::fixed<spacedim> >(this->n_dofs_*this->n_components_));
    }

    //views_cache_.initialize(*this->fv_, _fe);
}


template<unsigned int spacedim>
template<unsigned int DIM>
std::shared_ptr<FEInternalData> PatchFEValues<spacedim>::DimPatchFEValues::init_fe_data(const FiniteElement<DIM> &fe, const Quadrature &q)
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


template<unsigned int spacedim>
void PatchFEValues<spacedim>::DimPatchFEValues::reinit(PatchElementsList patch_elements)
{
    if (dim_ == 0) return; // Temporary skip, remove if PatchFEValues objects will be merge to common object of GenericAssembly

    element_patch_map_.clear();
    if (object_type_ == ElementFE)
        used_size_ = patch_elements.size();
    else
        used_size_ = patch_elements.size() * (this->dim_+1);
    ASSERT_LE(used_size_, max_size());

    unsigned int i=0;
    for (auto it=patch_elements.begin(); it!=patch_elements.end(); ++it, ++i) {
        if (object_type_ == ElementFE) {
            patch_data_idx_ = i;
            element_patch_map_[it->second] = i;
            element_data_[i].elm_values_->reinit(it->first);
            //this->fill_data(*element_data_[i].elm_values_, *this->fe_data_);
        } else {
            element_patch_map_[it->second] = i * (this->dim_+1);
            for (unsigned int sid=0; sid<this->dim_+1; ++sid) {
                patch_data_idx_ = i * (this->dim_+1) + sid;
                element_data_[patch_data_idx_].elm_values_->reinit( *it->first.side(sid) );
                //this->fill_data(*element_data_[patch_data_idx_].elm_values_, *this->side_fe_data_[sid]);

            }
        }
    }
}



// explicit instantiation
template void PatchFEValues<3>::initialize<0>(Quadrature&, FiniteElement<0>&, UpdateFlags);
template void PatchFEValues<3>::initialize<1>(Quadrature&, FiniteElement<1>&, UpdateFlags);
template void PatchFEValues<3>::initialize<2>(Quadrature&, FiniteElement<2>&, UpdateFlags);
template void PatchFEValues<3>::initialize<3>(Quadrature&, FiniteElement<3>&, UpdateFlags);

template class PatchFEValues<3>;
