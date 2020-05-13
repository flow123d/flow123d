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
#include "fem/fe_system.hh"



using namespace std;





template<class FS> double Dof::evaluate(const FS &function_space,
                                        unsigned int basis_idx) const
{
    // Check that FS is derived from FunctionSpace.
    static_assert(std::is_base_of<FunctionSpace, FS>::value, "FS must be derived from FunctionSpace.");
    
    // We cannot evaluate dof on dim-dimensional n-face if the function space lies on lower-dimensional n-face.
    ASSERT(function_space.space_dim()+1 == coords.size());
    
    switch (type)
    {
    case Value:
    {
        // evaluate basis function and return the linear combination of components
        arma::vec vec_value(function_space.n_components());
        
        // drop off the 0-th barycentric coordinate
        // in case of space_dim=0 subvec does not work
        arma::vec f_coords(function_space.space_dim());
        if (function_space.space_dim() > 0)
            f_coords = coords.subvec(1,coords.size()-1);
        for (unsigned int c=0; c<function_space.n_components(); c++)
            vec_value[c] = function_space.basis_value(basis_idx, f_coords, c);
        return dot(coefs, vec_value);
        break;
    }
        
    default:
        OLD_ASSERT(false, "Dof evaluation not implemented for this type.");
    }
    return 0;
}







template<unsigned int dim>
FiniteElement<dim>::FiniteElement()
    : function_space_(nullptr)
{
    init();
}

template<unsigned int dim>
void FiniteElement<dim>::init(bool primitive, FEType type)
{
    dofs_.clear();
    is_primitive_ = primitive;
    type_ = type;
}


template<unsigned int dim>
void FiniteElement<dim>::setup_components()
{
  component_indices_.resize(dofs_.size(), 0);
  nonzero_components_.resize(dofs_.size(), { true });
}


template<unsigned int dim> inline
void FiniteElement<dim>::compute_node_matrix()
{
    arma::mat M(dofs_.size(), dofs_.size());

    for (unsigned int i = 0; i < dofs_.size(); i++)
        for (unsigned int j = 0; j < dofs_.size(); j++) {
            M(j, i) = dofs_[i].evaluate(*function_space_, j);

        }
    node_matrix = arma::inv(M);
}


template<unsigned int dim>
double FiniteElement<dim>::shape_value(const unsigned int i, 
                                       const arma::vec::fixed<dim> &p,
                                       const unsigned int comp) const
{
    ASSERT_DBG( comp < n_components() );
	ASSERT_DBG( i < dofs_.size()).error("Index of basis function is out of range.");
    
    double value = 0;
    for (unsigned int j=0; j<function_space_->dim(); j++)
        value += function_space_->basis_value(j, p, comp) * node_matrix(i,j);

    return value;
}

template<unsigned int dim>
arma::vec::fixed<dim> FiniteElement<dim>::shape_grad(const unsigned int i,
                                                     const arma::vec::fixed<dim> &p,
                                                     const unsigned int comp) const
{
    ASSERT_DBG( comp < n_components() );
	ASSERT_DBG( i < dofs_.size()).error("Index of basis function is out of range.");
    
    // uword is a typedef for an unsigned integer type; it is used for matrix indices as well as all internal counters and loops
    // sword is a typedef for a signed integer type
    arma::vec grad((arma::uword)dim);
    grad.zeros();
    for (unsigned int j=0; j<function_space_->dim(); j++)
        grad += function_space_->basis_grad(j, p, comp) * node_matrix(i,j);
    
    return grad;
}


template<unsigned int dim> inline
UpdateFlags FiniteElement<dim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    switch (type_)
    {
        case FEScalar:   
        case FEVector:
        case FETensor:
            if (flags & update_gradients)
                f |= update_inverse_jacobians;
            break;
        case FEVectorContravariant:
            if (flags & update_values)
                f |= update_jacobians;
            if (flags & update_gradients)
                f |= update_jacobians | update_inverse_jacobians;
            break;
        case FEVectorPiola:
            if (flags & update_values)
                f |= update_jacobians | update_volume_elements;
            if (flags & update_gradients)
                f |= update_jacobians | update_inverse_jacobians | update_volume_elements;
            break;
        default:;
    }

    return f;
}


template<unsigned int dim>
unsigned int FiniteElement<dim>::n_space_components(unsigned int spacedim)
{
    switch (type_) {
        case FEScalar:
            return 1;
            break;
        case FEVector:
        case FEVectorContravariant:
        case FEVectorPiola:
            return spacedim;
            break;
        case FETensor:
            return spacedim*spacedim;
            break;
        case FEMixedSystem:
            const FESystem<dim> *fe_sys = dynamic_cast<const FESystem<dim>*>(this);
            ASSERT_DBG(fe_sys != nullptr).error("Mixed system must be represented by FESystem.");
            return  fe_sys->get_scalar_components().size()
                   +fe_sys->get_vector_components().size()*spacedim
                   +fe_sys->get_tensor_components().size()*spacedim*spacedim;
            break;
    }

    // should be never reached
    ASSERT(0).error("Unknown type of FiniteElement.");
    return 0;
}


template<unsigned int dim>
vector< arma::vec::fixed<dim+1> > FiniteElement<dim>::dof_points() const {
    std::vector<arma::vec::fixed<dim+1>> points(20);
    points.resize(0);
    for(auto dof : this->dofs_)
        points.push_back(dof.coords);
    return points;
}





template class FiniteElement<0>;
template class FiniteElement<1>;
template class FiniteElement<2>;
template class FiniteElement<3>;


