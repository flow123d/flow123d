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

#include "fem/mapping.hh"
#include "quadrature/quadrature.hh"
#include <quadrature/qxfem.hh>
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"
#include "fem/fe_system.hh"



using namespace arma;
using namespace std;








FEInternalData::FEInternalData(unsigned int np, unsigned int nd)
    : n_points(np),
      n_dofs(nd)
{
    ref_shape_values.resize(np, vector<arma::vec>(nd));
    ref_shape_grads.resize(np, vector<arma::mat>(nd));
}



template<unsigned int dim, unsigned int spacedim>
void FEValuesData<dim,spacedim>::allocate(unsigned int size, UpdateFlags flags, unsigned int n_comp)
{
    update_flags = flags;

    // resize the arrays of computed quantities
    if (update_flags & update_jacobians)
        jacobians.resize(size);

    if (update_flags & update_volume_elements)
        determinants.resize(size);

    if ((update_flags & update_JxW_values) |
        (update_flags & update_side_JxW_values))
        JxW_values.resize(size);

    if (update_flags & update_inverse_jacobians)
        inverse_jacobians.resize(size);

    if (update_flags & update_values)
    {
        shape_values.clear();
        shape_values.resize(size, vector<double>(n_comp));
    }

    if (update_flags & update_gradients)
    {
        shape_gradients.clear();
        shape_gradients.resize(size, vector<arma::vec::fixed<spacedim> >(n_comp));
    }
    
    if (update_flags & update_divergence){
        shape_divergence.resize(size);
    }

    if (update_flags & update_quadrature_points)
        points.resize(size);

    if (update_flags & update_normal_vectors)
        normal_vectors.resize(size);
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
      FESystem<dim> *fe = dynamic_cast<FESystem<dim>*>(fv.get_fe());
      OLD_ASSERT(fe != nullptr, "Mixed system must be represented by FESystem.");
      std::vector<unsigned int> sc = fe->get_scalar_components();
      std::vector<unsigned int> vc = fe->get_vector_components();
      std::vector<unsigned int> tc = fe->get_tensor_components();
      for (auto si : sc)
        scalars.push_back(FEValuesViews::Scalar<dim,spacedim>(fv,si));
      for (auto vi : vc)
        vectors.push_back(FEValuesViews::Vector<dim,spacedim>(fv,vi));
      for (auto ti : tc)
        tensors.push_back(FEValuesViews::Tensor<dim,spacedim>(fv,ti));
      break;
  }
}



template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::FEValuesBase()
: mapping(NULL), quadrature(NULL), fe(NULL), mapping_data(NULL), fe_data(NULL)
{
}



template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::~FEValuesBase() {
    if (mapping_data) delete mapping_data;
    if (fe_data) delete fe_data;
}



template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::allocate(Mapping<dim,spacedim> & _mapping,
        Quadrature<dim> & _quadrature,
        FiniteElement<dim> & _fe,
        unsigned int ndofs,
        UpdateFlags _flags)
{
    // For FEVector and FETensor check number of components.
    // This cannot be done in FiniteElement since it does not know spacedim.
    if (_fe.type_ == FEVector)
        ASSERT_DBG(_fe.n_components() == spacedim).error("FEVector must have spacedim components.");
    else if (_fe.type_ == FETensor)
        ASSERT_DBG(_fe.n_components() == spacedim*spacedim).error("FETensor must have spacedim*spacedim components.");
    
    mapping = &_mapping;
    quadrature = &_quadrature;
    fe = &_fe;
    n_dofs_ = ndofs;

    switch (fe->type_) {
        case FEScalar:
            n_components_ = 1;
            break;
        case FEVector:
        case FEVectorContravariant:
        case FEVectorPiola:
            n_components_ = spacedim;
            break;
        case FETensor:
            n_components_ = spacedim*spacedim;
            break;
        case FEMixedSystem:
            n_components_ = fe->n_components();
            break;
    }
    // add flags required by the finite element or mapping
    data.allocate(quadrature->size(), update_each(_flags), ndofs*n_components_);
    
    views_cache_.initialize(*this);
}



template<unsigned int dim, unsigned int spacedim>
FEInternalData *FEValuesBase<dim,spacedim>::init_fe_data(const Quadrature<dim> *q)
{
    FEInternalData *data = new FEInternalData(q->size(), fe->n_dofs());

    arma::mat shape_values(fe->n_dofs(), fe->n_components());
    for (unsigned int i=0; i<q->size(); i++)
    {
        for (unsigned int j=0; j<fe->n_dofs(); j++)
        {
            for (unsigned int c=0; c<fe->n_components(); c++)
                shape_values(j,c) = fe->shape_value(j, q->point(i), c);
            
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
                grad.col(c) += fe->shape_grad(j, q->point(i), c);
            
            data->ref_shape_grads[i][j] = grad;
        }
    }
    
    return data;
}



template<unsigned int dim, unsigned int spacedim>
UpdateFlags FEValuesBase<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags | fe->update_each(flags);
    f |= mapping->update_each(f);
    return f;
}


template<unsigned int dim, unsigned int spacedim>
double FEValuesBase<dim,spacedim>::shape_value(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, quadrature->size());
  return data.shape_values[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_grad(const unsigned int function_no, const unsigned int point_no)
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, quadrature->size());
  return data.shape_gradients[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim>
double FEValuesBase<dim,spacedim>::shape_value_component(const unsigned int function_no, 
                                    const unsigned int point_no, 
                                    const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, quadrature->size());
  ASSERT_LT_DBG(comp, n_components_);
  return data.shape_values[point_no][function_no*n_components_+comp];
}


template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_grad_component(const unsigned int function_no,
                                                        const unsigned int point_no,
                                                        const unsigned int comp) const
{
  ASSERT_LT_DBG(function_no, n_dofs_);
  ASSERT_LT_DBG(point_no, quadrature->size());
  ASSERT_LT_DBG(comp, n_components_);
  return data.shape_gradients[point_no][function_no*n_components_+comp];
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_scalar_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEScalar);
    
    // shape values
    if (data.update_flags & update_values)
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                data.shape_values[i][j] = fe_data.ref_shape_values[i][j][0];

    // shape gradients
    if (data.update_flags & update_gradients)
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                data.shape_gradients[i][j] = trans(data.inverse_jacobians[i]) * fe_data.ref_shape_grads[i][j];
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVector);
    
    // shape values
    if (data.update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (data.update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(data.inverse_jacobians[i]) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_contravariant_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVectorContravariant);
    
    // shape values
    if (data.update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = data.jacobians[i] * fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (data.update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(data.inverse_jacobians[i]) * fe_data.ref_shape_grads[i][j] * trans(data.jacobians[i]);
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_vec_piola_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEVectorPiola);
    
    // shape values
    if (data.update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = data.jacobians[i]*fe_data.ref_shape_values[i][j]/data.determinants[i];
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_values[i][j*spacedim+c] = fv_vec(c);
            }
    }

    // shape gradients
    if (data.update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(data.inverse_jacobians[i]) * fe_data.ref_shape_grads[i][j] * trans(data.jacobians[i])
                        / data.determinants[i];
                for (unsigned int c=0; c<spacedim; c++)
                    data.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }   
    }
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_tensor_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FETensor);
    
    // shape values
    if (data.update_flags & update_values)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    data.shape_values[i][j*spacedim*spacedim+c] = fv_vec[c];
            }
    }

    // shape gradients
    if (data.update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(data.inverse_jacobians[i]) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    data.shape_gradients[i][j*spacedim*spacedim+c] = grads.col(c);
            }
    }
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_system_data(const FEInternalData &fe_data)
{
    ASSERT_DBG(fe->type_ == FEMixedSystem);
    
    // for mixed system we first fill data in sub-elements
    FESystem<dim> *fe_sys = dynamic_cast<FESystem<dim>*>(fe);
    ASSERT_DBG(fe_sys != nullptr).error("Mixed system must be represented by FESystem.");
    for (unsigned int f=0; f<fe_sys->fe().size(); f++)
    {
        // fill fe_values for base FE
        fe_values_vec[f]->data.jacobians = data.jacobians;
        fe_values_vec[f]->data.inverse_jacobians = data.inverse_jacobians;
        fe_values_vec[f]->data.determinants = data.determinants;
        fe_values_vec[f]->fill_data(*fe_values_vec[f]->fe_data);
    }
    
    // shape values
    if (data.update_flags & update_values)
    {
        arma::vec fv_vec;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys->fe().size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys->fe()[f]->n_dofs(); n++)
                    for (unsigned int c=0; c<fe_sys->fe()[f]->n_components(); c++)
                        data.shape_values[i][shape_offset+fe_sys->function_space_->n_components()*n+comp_offset+c] = fe_values_vec[f]->data.shape_values[i][n*fe_sys->fe()[f]->n_components()+c];
            
            comp_offset += fe_sys->fe()[f]->n_components();
            shape_offset += fe_sys->fe()[f]->n_dofs()*fe_sys->function_space_->n_components();
        }
    }

    // shape gradients
    if (data.update_flags & update_gradients)
    {
        arma::mat grads;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_sys->fe().size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_sys->fe()[f]->n_dofs(); n++)
                    for (unsigned int c=0; c<fe_sys->fe()[f]->n_components(); c++)
                        data.shape_gradients[i][shape_offset+fe_sys->function_space_->n_components()*n+comp_offset+c] = fe_values_vec[f]->data.shape_gradients[i][n*fe_sys->fe()[f]->n_components()+c];
            
            comp_offset += fe_sys->fe()[f]->n_components();
            shape_offset += fe_sys->fe()[f]->n_dofs()*fe_sys->function_space_->n_components();
        }
    }
    
}


template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::fill_data(const FEInternalData &fe_data)
{
    switch (fe->type_) {
        case FEScalar:
            fill_scalar_data(fe_data);
            break;
        case FEVector:
            fill_vec_data(fe_data);
            break;
        case FEVectorContravariant:
            fill_vec_contravariant_data(fe_data);
            break;
        case FEVectorPiola:
            fill_vec_piola_data(fe_data);
            break;
        case FETensor:
            fill_tensor_data(fe_data);
            break;
        case FEMixedSystem:
            fill_system_data(fe_data);
            break;
        default:
            ASSERT(false).error("Not implemented.");
    }
}








template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::FEValues(Mapping<dim,spacedim> &_mapping,
         Quadrature<dim> &_quadrature,
         FiniteElement<dim> &_fe,
         UpdateFlags _flags)
:FEValuesBase<dim, spacedim>()
{
    this->allocate(_mapping, _quadrature, _fe, _fe.n_dofs(), _flags);

    //UpdateFlags map_update_flags = this->data.update_flags;
    
    if(typeid(_quadrature) == typeid(QXFEM<dim,spacedim>)){
//         DBGCOUT("FEValues: quad XFEM - dismiss update_quadrature_points flag.\n");
        QXFEM<dim,spacedim>* q = static_cast<QXFEM<dim,spacedim>*>(&_quadrature);
        
//         map_update_flags = map_update_flags & (~update_quadrature_points);
        this->data.update_flags = this->data.update_flags & (~update_quadrature_points);
        //TODO: think of way to avoid copying the whole vector
        this->data.points = q->get_real_points();
    }
    
    // precompute the maping data and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    this->fe_data = this->init_fe_data(this->quadrature);
    
    // In case of mixed system allocate data for sub-elements.
    if (this->fe->type_ == FEMixedSystem)
    {
        FESystem<dim> *fe = dynamic_cast<FESystem<dim>*>(this->fe);
        ASSERT_DBG(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        for (auto fe_sub : fe->fe())
            this->fe_values_vec.push_back(make_shared<FEValues<dim,spacedim> >(_mapping, _quadrature, *fe_sub, _flags));
    }
}



template<unsigned int dim,unsigned int spacedim>
void FEValues<dim,spacedim>::reinit(ElementAccessor<3> & cell)
{
	OLD_ASSERT_EQUAL( dim, cell->dim() );
    this->data.present_cell = &cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_values(cell,
                            *this->quadrature,
                            *this->mapping_data,
                            this->data);

    this->fill_data(*this->fe_data);
}









template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::FESideValues(Mapping<dim,spacedim> & _mapping,
                                 Quadrature<dim-1> & _sub_quadrature,
                                 FiniteElement<dim> & _fe,
                                 const UpdateFlags _flags)
:FEValuesBase<dim,spacedim>()
{
    sub_quadrature = &_sub_quadrature;
    Quadrature<dim> *q = new Quadrature<dim>(_sub_quadrature.size());
    this->allocate(_mapping, *q, _fe, _fe.n_dofs(), _flags);

    for (unsigned int sid = 0; sid < RefElement<dim>::n_sides; sid++)
    {
    	for (unsigned int pid = 0; pid < RefElement<dim>::n_side_permutations; pid++)
    	{
    		// transform the side quadrature points to the cell quadrature points
            side_quadrature[sid][pid] = Quadrature<dim>(_sub_quadrature, sid, pid);
    		side_mapping_data[sid][pid] = this->mapping->initialize(side_quadrature[sid][pid], this->data.update_flags);
    		side_fe_data[sid][pid] = this->init_fe_data(&side_quadrature[sid][pid]);
    	}
    }
    
    // In case of mixed system allocate data for sub-elements.
    if (this->fe->type_ == FEMixedSystem)
    {
        FESystem<dim> *fe = dynamic_cast<FESystem<dim>*>(this->fe);
        ASSERT_DBG(fe != nullptr).error("Mixed system must be represented by FESystem.");
        
        for (auto fe_sub : fe->fe())
            this->fe_values_vec.push_back(make_shared<FESideValues<dim,spacedim> >(_mapping, _sub_quadrature, *fe_sub, _flags));
    }
}



template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::~FESideValues()
{
	for (unsigned int sid=0; sid<RefElement<dim>::n_sides; sid++)
	{
		for (unsigned int pid=0; pid<RefElement<dim>::n_side_permutations; pid++)
		{
			delete side_mapping_data[sid][pid];
			delete side_fe_data[sid][pid];
		}
	}

    // Since quadrature is an auxiliary internal variable allocated
    // by the constructor, it must be destroyed here.
    delete this->quadrature;
}


template<unsigned int dim,unsigned int spacedim>
void FESideValues<dim,spacedim>::reinit(ElementAccessor<3> & cell,
		unsigned int sid)
{
    ASSERT_LT_DBG( sid, cell->n_sides());
    ASSERT_EQ_DBG(dim, cell->dim());
    this->data.present_cell = &cell;

    unsigned int pid = cell->permutation_idx(sid);
    ASSERT_LT_DBG(pid, RefElement<dim>::n_side_permutations);
    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_side_values(cell,
                                 sid,
                                 side_quadrature[sid][pid],
                                 *side_mapping_data[sid][pid],
                                 this->data);

    // calculation of finite element data
    this->fill_data(*side_fe_data[sid][pid]);
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

