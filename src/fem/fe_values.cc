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







template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::FEInternalData::FEInternalData(unsigned int np, unsigned int nd)
    : n_points(np),
      n_dofs(nd)
{
    ref_shape_values.resize(np, vector<arma::vec>(nd));
    ref_shape_grads.resize(np, vector<arma::mat>(nd));
}


template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::FEInternalData::FEInternalData(const FEInternalData &fe_system_data,
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



template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::ViewsCache::initialize(FEValuesBase<dim,spacedim> &fv)
{
  scalars.clear();
  vectors.clear();
  tensors.clear();
  switch (fv.get_fe()->type_) {
    case FEType::FEScalar:
      scalars.push_back(FEValuesViews::Scalar<dim,spacedim>(fv, 0));
      break;
    case FEType::FEVector:
    case FEType::FEVectorContravariant:
    case FEType::FEVectorPiola:
      vectors.push_back(FEValuesViews::Vector<dim,spacedim>(fv, 0));
      break;
    case FEType::FETensor:
      tensors.push_back(FEValuesViews::Tensor<dim,spacedim>(fv, 0));
      break;
    case FEType::FEMixedSystem:
      FESystem<dim> *fe_sys = dynamic_cast<FESystem<dim>*>(fv.get_fe());
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
              scalars.push_back(FEValuesViews::Scalar<dim,spacedim>(fv,comp_offset));
              break;
          case FEType::FEVector:
          case FEType::FEVectorContravariant:
          case FEType::FEVectorPiola:
              vectors.push_back(FEValuesViews::Vector<dim,spacedim>(fv,comp_offset));
              break;
          case FEType::FETensor:
              tensors.push_back(FEValuesViews::Tensor<dim,spacedim>(fv,comp_offset));
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



template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::FEValuesBase()
: n_points_(0), fe(nullptr), elm_values(nullptr)
{
}



template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::~FEValuesBase() {
    if (elm_values != nullptr) delete elm_values;
}



template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::allocate(
        unsigned int n_points,
        FiniteElement<dim> & _fe,
        UpdateFlags _flags)
{
    // For FEVector and FETensor check number of components.
    // This cannot be done in FiniteElement since it does not know spacedim.
    if (_fe.type_ == FEVector) {
        ASSERT_DBG(_fe.n_components() == spacedim).error("FEVector must have spacedim components.");
    } else if (_fe.type_ == FETensor) {
        ASSERT_DBG(_fe.n_components() == spacedim*spacedim).error("FETensor must have spacedim*spacedim components.");
    }
    
    n_points_ = n_points;
    fe = &_fe;
    n_components_ = fe->n_space_components(spacedim);
    
    // add flags required by the finite element or mapping
    update_flags = update_each(_flags);
    if (update_flags & update_values)
        shape_values.resize(n_points_, vector<double>(fe->n_dofs()*n_components_));

    if (update_flags & update_gradients)
        shape_gradients.resize(n_points_, vector<arma::vec::fixed<spacedim> >(fe->n_dofs()*n_components_));
    
    views_cache_.initialize(*this);
}



template<unsigned int dim, unsigned int spacedim>
typename FEValuesBase<dim,spacedim>::FEInternalData *FEValuesBase<dim,spacedim>::init_fe_data(const Quadrature *q)
{
    ASSERT_DBG( q->dim() == dim );
    FEInternalData *data = new FEInternalData(q->size(), fe->n_dofs());

    //DebugOut() << "q size: " << q->size() << "\n";
    arma::mat shape_values(fe->n_dofs(), fe->n_components());
    for (unsigned int i=0; i<q->size(); i++)
    {
        for (unsigned int j=0; j<fe->n_dofs(); j++)
        {
            for (unsigned int c=0; c<fe->n_components(); c++)
                shape_values(j,c) = fe->shape_value(j, q->point<dim>(i), c);
            
            data->ref_shape_values[i][j] = trans(shape_values.row(j));
        }
    }

    arma::mat grad(dim, fe->n_components());
    for (unsigned int i=0; i<q->size(); i++)
    {
        for (unsigned int j=0; j<fe->n_dofs(); j++)
        {
            grad.zeros();
            for (unsigned int c=0; c<fe->n_components(); c++)
                grad.col(c) += fe->shape_grad(j, q->point<dim>(i), c);
            
            data->ref_shape_grads[i][j] = grad;
        }
    }
    
    return data;
}



template<unsigned int dim, unsigned int spacedim>
UpdateFlags FEValuesBase<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags | fe->update_each(flags);
    f |= MappingP1<dim,spacedim>::update_each(f);
    return f;
}


template<unsigned int dim, unsigned int spacedim>
double FEValuesBase<dim,spacedim>::shape_value(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, fe->n_dofs());
  ASSERT_LT_DBG(point_no, n_points_);
  return shape_values[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_grad(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, fe->n_dofs());
  ASSERT_LT_DBG(point_no, n_points_);
  return shape_gradients[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim>
double FEValuesBase<dim,spacedim>::shape_value_component(const unsigned int function_no, 
                                    const unsigned int point_no, 
                                    const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, fe->n_dofs());
  ASSERT_LT_DBG(point_no, n_points_);
  ASSERT_LT_DBG(comp, n_components_);
  return shape_values[point_no][function_no*n_components_+comp];
}


template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_grad_component(const unsigned int function_no,
                                                        const unsigned int point_no,
                                                        const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, fe->n_dofs());
  ASSERT_LT_DBG(point_no, n_points_);
  ASSERT_LT_DBG(comp, n_components_);
  return shape_gradients[point_no][function_no*n_components_+comp];
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_scalar_data(const ElementValuesBase<dim,spacedim> &elm_values, const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEScalar);
    
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


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_data(const ElementValuesBase<dim,spacedim> &elm_values,
                                               const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVector);
    
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


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_contravariant_data(const ElementValuesBase<dim,spacedim> &elm_values,
                                                             const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVectorContravariant);
    
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


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_piola_data(const ElementValuesBase<dim,spacedim> &elm_values,
                                                     const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVectorPiola);
    
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


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_tensor_data(const ElementValuesBase<dim,spacedim> &elm_values,
                                                  const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FETensor);
    
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


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_system_data(const ElementValuesBase<dim,spacedim> &elm_values, const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEMixedSystem);
    
    // for mixed system we first fill data in sub-elements
    FESystem<dim> *fe_sys = dynamic_cast<FESystem<dim>*>(fe);
    ASSERT_DBG(fe_sys != nullptr).error("Mixed system must be represented by FESystem.");
    FESystemFunctionSpace *fs = dynamic_cast<FESystemFunctionSpace*>(fe_sys->function_space_.get());
    unsigned int comp_offset = 0;
    for (unsigned int f=0; f<fe_sys->fe().size(); f++)
    {
        // fill fe_values for base FE
        unsigned int n_comp = fe_sys->fe()[f]->n_components();
        FEInternalData vec_fe_data(fe_data, fe_sys->fe_dofs(f), comp_offset, n_comp);
        fe_values_vec[f]->fill_data(elm_values, vec_fe_data);
        
        comp_offset += n_comp;
    }
    
    unsigned int n_space_components = fe->n_space_components(spacedim);
    
    // shape values
    if (update_flags & update_values)
    {
        arma::vec fv_vec;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys->fe().size(); f++)
        {
            unsigned int n_sub_space_components = fe_sys->fe()[f]->n_space_components(spacedim);
            
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys->fe()[f]->n_dofs(); n++)
                    for (unsigned int c=0; c<n_sub_space_components; c++)
                        shape_values[i][shape_offset+n_space_components*n+comp_offset+c] = fe_values_vec[f]->shape_values[i][n*n_sub_space_components+c];
            
            comp_offset += n_sub_space_components;
            shape_offset += fe_sys->fe()[f]->n_dofs()*n_space_components;
        }
    }

    // shape gradients
    if (update_flags & update_gradients)
    {
        arma::mat grads;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys->fe().size(); f++)
        {
            unsigned int n_sub_space_components = fe_sys->fe()[f]->n_space_components(spacedim);
            
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys->fe()[f]->n_dofs(); n++)
                    for (unsigned int c=0; c<n_sub_space_components; c++)
                        shape_gradients[i][shape_offset+n_space_components*n+comp_offset+c] = fe_values_vec[f]->shape_gradients[i][n*n_sub_space_components+c];
            
            comp_offset += n_sub_space_components;
            shape_offset += fe_sys->fe()[f]->n_dofs()*n_space_components;
        }
    }
    
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_data(const ElementValuesBase<dim,spacedim> &elm_values, const FEInternalData &fe_data)
{
    switch (fe->type_) {
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








template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::FEValues(
         Quadrature &q,
         FiniteElement<dim> &_fe,
         UpdateFlags _flags)
: FEValuesBase<dim, spacedim>(),
  fe_data(nullptr)
{
    if (dim == 0) return; // avoid unnecessary allocation of dummy 0 dimensional objects
    ASSERT_DBG( q.dim() == dim );
    this->allocate( q.size(), _fe, _flags);
    this->elm_values = new ElementValues<dim,spacedim>(q, this->update_flags);
    
    // precompute finite element data
    fe_data = this->init_fe_data(&q);
    
    // In case of mixed system allocate data for sub-elements.
    if (this->fe->type_ == FEMixedSystem)
    {
        FESystem<dim> *fe = dynamic_cast<FESystem<dim>*>(this->fe);
        ASSERT_DBG(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        for (auto fe_sub : fe->fe())
            this->fe_values_vec.push_back(make_shared<FEValues<dim,spacedim> >(q, *fe_sub, this->update_flags));
    }
}


template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::~FEValues()
{
    if (fe_data) delete fe_data;
}



template<unsigned int dim,unsigned int spacedim>
void FEValues<dim,spacedim>::reinit(const ElementAccessor<spacedim> &cell)
{
	OLD_ASSERT_EQUAL( dim, cell->dim() );
    
    if (!this->elm_values->cell().is_valid() ||
        this->elm_values->cell().idx() != cell.idx())
    {
        ((ElementValues<dim,spacedim> *)this->elm_values)->reinit(cell);
    }
    
    this->fill_data(*this->elm_values, *fe_data);
}



MixedPtr<FEValues> mixed_fe_values(
        QGauss::array &quadrature,
        MixedPtr<FiniteElement> fe,
        UpdateFlags flags)
{
    return MixedPtr<FEValues>(
      std::make_shared<FEValues<0>>(quadrature[0], *fe.get<0>(), flags),
      std::make_shared<FEValues<1>>(quadrature[1], *fe.get<1>(), flags),
      std::make_shared<FEValues<2>>(quadrature[2], *fe.get<2>(), flags),
      std::make_shared<FEValues<3>>(quadrature[3], *fe.get<3>(), flags)
      );
}






template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::FESideValues(
                                 Quadrature & _sub_quadrature,
                                 FiniteElement<dim> & _fe,
                                 const UpdateFlags _flags)
: FEValuesBase<dim,spacedim>(),
  side_idx_(-1)
{
    ASSERT_DBG( _sub_quadrature.dim() + 1 == dim );
    this->allocate( _sub_quadrature.size(), _fe, _flags);
    this->elm_values = new ElemSideValues<dim,spacedim>( _sub_quadrature, this->update_flags );
    
    for (unsigned int sid = 0; sid < RefElement<dim>::n_sides; sid++)
    {
    	for (unsigned int pid = 0; pid < RefElement<dim>::n_side_permutations; pid++)
    	{
    		// transform the side quadrature points to the cell quadrature points
    		side_fe_data[sid][pid] = this->init_fe_data(&((ElemSideValues<dim,spacedim> *)this->elm_values)->quadrature(sid,pid));
    	}
    }
    

    // In case of mixed system allocate data for sub-elements.
    if (this->fe->type_ == FEMixedSystem)
    {
        FESystem<dim> *fe = dynamic_cast<FESystem<dim>*>(this->fe);
        ASSERT_DBG(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        for (auto fe_sub : fe->fe())
            this->fe_values_vec.push_back(make_shared<FESideValues<dim,spacedim> >(_sub_quadrature, *fe_sub, this->update_flags));
    }
}



template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::~FESideValues()
{
    for (unsigned int sid=0; sid<RefElement<dim>::n_sides; sid++)
        for (unsigned int pid=0; pid<RefElement<dim>::n_side_permutations; pid++)
            delete side_fe_data[sid][pid];
}


template<unsigned int dim,unsigned int spacedim>
void FESideValues<dim,spacedim>::reinit(const ElementAccessor<spacedim> &cell, unsigned int sid)
{
    ASSERT_LT_DBG( sid, cell->n_sides() );
    ASSERT_EQ_DBG( dim, cell->dim() );
    
    if (!this->elm_values->cell().is_valid() || 
        this->elm_values->cell().idx() != cell.idx() ||
        side_idx_ != sid)
    {
        ((ElemSideValues<dim,spacedim> *)this->elm_values)->reinit(cell, sid);
    }
    
    side_idx_ = sid;
    
    // calculation of finite element data
    this->fill_data(*this->elm_values, *side_fe_data[sid][cell->permutation_idx(sid)]);
}







template class FEValuesBase<0,3>;
template class FEValuesBase<1,3>;
template class FEValuesBase<2,3>;
template class FEValuesBase<3,3>;

template class FEValues<0,3>;
template class FEValues<1,3>;
template class FEValues<2,3>;
template class FEValues<3,3>;

template class FESideValues<1,3>;
template class FESideValues<2,3>;
template class FESideValues<3,3>;

