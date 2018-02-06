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
 * @file    finite_element.cc
 * @brief   Abstract class for description of finite elements.
 * @author  Jan Stebel
 */

#include "system/system.hh"
#include "quadrature/quadrature.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"



using namespace std;





template<class FS> const double Dof::evaluate(const FS &function_space,
                                              unsigned int basis_idx) const
{
    // We cannot evaluate dof on dim-dimensional n-face if the function space lies on lower-dimensional n-face.
    ASSERT(function_space.space_dim()+1 == coords.size());
    
    switch (type)
    {
    case Value:
    {
        // evaluate basis function and return the linear combination of components
        arma::vec vec_value(function_space.n_components());
        for (unsigned int c=0; c<function_space.n_components(); c++)
            vec_value[c] = function_space.basis_value(basis_idx, coords.subvec(0,coords.size()-2), c);
        return dot(coefs, vec_value);
        break;
    }
        
    default:
        OLD_ASSERT(false, "Dof evaluation not implemented for this type.");
    }
    return 0;
}



template<unsigned int dim, unsigned int spacedim>
FiniteElement<dim,spacedim>::FiniteElement()
    : function_space_(nullptr)
{
    init();
}

template<unsigned int dim, unsigned int spacedim>
void FiniteElement<dim,spacedim>::init(unsigned int n_components, bool primitive, FEType type)
{
    dofs_.clear();
    for (unsigned int i = 0; i <= dim; i++)
    {
        number_of_single_dofs[i] = 0;
        number_of_pairs[i] = 0;
        number_of_triples[i] = 0;
        number_of_sextuples[i] = 0;
    }
    
    is_primitive_ = primitive;
    n_components_ = n_components;
    type_ = type;
}


template<unsigned int dim, unsigned int spacedim>
void FiniteElement<dim,spacedim>::setup_components()
{
  component_indices_.resize(dofs_.size(), 0);
  nonzero_components_.resize(dofs_.size(), { true });
}


template<unsigned int dim, unsigned int spacedim> inline
void FiniteElement<dim,spacedim>::compute_node_matrix()
{
    arma::mat M(dofs_.size(), dofs_.size());

    for (unsigned int i = 0; i < dofs_.size(); i++)
        for (unsigned int j = 0; j < dofs_.size(); j++) {
            M(j, i) = dofs_[i].evaluate(*function_space_, j);

        }
    node_matrix = arma::inv(M);
}

template<unsigned int dim, unsigned int spacedim>
FEInternalData *FiniteElement<dim,spacedim>::initialize(const Quadrature<dim> &q)
{
    FEInternalData *data = new FEInternalData;

    arma::vec values(dofs_.size());
    data->basis_values.resize(q.size());
    for (unsigned int i=0; i<q.size(); i++)
    {
        for (unsigned int j=0; j<dofs_.size(); j++)
            values[j] = basis_value(j, q.point(i));
        data->basis_values[i] = node_matrix * values;
    }

    arma::mat grads(dofs_.size(), dim);
    data->basis_grads.resize(q.size());
    for (unsigned int i=0; i<q.size(); i++)
    {
        for (unsigned int j=0; j<dofs_.size(); j++)
            grads.row(j) = arma::trans(basis_grad(j, q.point(i)));
        data->basis_grads[i] = node_matrix * grads;
    }

    return data;
}


template<unsigned int dim, unsigned int spacedim>
double FiniteElement<dim,spacedim>::basis_value(const unsigned int i, 
                                       const arma::vec::fixed<dim> &p,
                                       const unsigned int comp) const
{
    ASSERT_DBG( comp < n_components_ );
	ASSERT_DBG( i < dofs_.size()).error("Index of basis function is out of range.");
    return this->function_space_->basis_value(i, p, comp);
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FiniteElement<dim,spacedim>::basis_grad(const unsigned int i,
                                                     const arma::vec::fixed<dim> &p,
                                                     const unsigned int comp) const
{
    ASSERT_DBG( comp < n_components_ );
	ASSERT_DBG( i < dofs_.size()).error("Index of basis function is out of range.");
    return this->function_space_->basis_grad(i, p, comp);
}


template<unsigned int dim, unsigned int spacedim> inline
UpdateFlags FiniteElement<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_gradients)
        f |= update_inverse_jacobians;

    return f;
}

template<unsigned int dim, unsigned int spacedim> inline
void FiniteElement<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &q,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    // shape values
    if (fv_data.update_flags & update_values)
    {
        for (unsigned int i = 0; i < q.size(); i++)
            for (unsigned int c = 0; c < n_dofs(); c++)
                fv_data.shape_values[i][c] = data.basis_values[i][c];
    }

    // shape gradients
    if (fv_data.update_flags & update_gradients)
    {
        for (unsigned int i = 0; i < q.size(); i++)
        {
            arma::mat grads = trans(data.basis_grads[i] * fv_data.inverse_jacobians[i]);
            for (unsigned int c = 0; c < n_dofs(); c++)
                fv_data.shape_gradients[i][c] = grads.col(c);
        }
    }
}



template<unsigned int dim, unsigned int spacedim>
FiniteElement<dim,spacedim>::~FiniteElement()
{
    if (function_space_ != nullptr) delete function_space_;
}


template class FiniteElement<0,3>;
template class FiniteElement<1,3>;
template class FiniteElement<2,3>;
template class FiniteElement<3,3>;


