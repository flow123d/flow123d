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
#include "fem/fe_values_views.hh"
#include "fem/fe_values.hh"
#include "fem/finite_element.hh"
#include "quadrature/quadrature.hh"

using namespace FEValuesViews;

template<unsigned int spacedim>
double FEValuesViews::Scalar<spacedim>::value(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  return fe_values_.shape_value_component(function_no, point_no, component_);
}

template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesViews::Scalar<spacedim>::grad(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  return fe_values_.shape_grad_component(function_no, point_no, component_);
}

template<unsigned int spacedim>
const FEValues<spacedim> &FEValuesViews::Scalar<spacedim>::base() const
{ return fe_values_; }
  



template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesViews::Vector<spacedim>::value(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::vec::fixed<spacedim> v;
  for (unsigned int c=0; c<spacedim; ++c)
    v(c) = fe_values_.shape_value_component(function_no, point_no, first_vector_component_+c);
  return v;
}

template<unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Vector<spacedim>::grad(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::mat::fixed<spacedim,spacedim> g;
  for (unsigned int c=0; c<spacedim; ++c)
    g.col(c) = fe_values_.shape_grad_component(function_no, point_no, first_vector_component_+c);
  return g.t();
}

template<unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Vector<spacedim>::sym_grad(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::mat::fixed<spacedim,spacedim> g = grad(function_no, point_no);
  return 0.5*(g+trans(g));
}

template<unsigned int spacedim>
double FEValuesViews::Vector<spacedim>::divergence(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  double div = 0;
  for (unsigned int c=0; c<spacedim; ++c)
    div += fe_values_.shape_grad_component(function_no, point_no, first_vector_component_+c)(first_vector_component_+c);
  return div;
}

template<unsigned int spacedim>
const FEValues<spacedim> &FEValuesViews::Vector<spacedim>::base() const
{ return fe_values_; }



template<unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Tensor<spacedim>::value(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::mat::fixed<spacedim,spacedim> v;
  for (unsigned int c=0; c<spacedim*spacedim; ++c)
    v(c/spacedim,c%spacedim) = fe_values_.shape_value_component(function_no, point_no, first_tensor_component_+c);
  return v;
}

template<unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Tensor<spacedim>::derivative(
    unsigned int variable_no,
    unsigned int function_no, 
    unsigned int point_no) const
{
  ASSERT_LT_DBG( variable_no, spacedim );
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::mat::fixed<spacedim,spacedim> g;
  for (unsigned int c=0; c<spacedim*spacedim; ++c)
    g(c/spacedim,c%spacedim) = fe_values_.shape_grad_component(function_no, point_no, first_tensor_component_+c)[first_tensor_component_+variable_no];
  return g;
}

template<unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesViews::Tensor<spacedim>::divergence(unsigned int function_no, unsigned int point_no) const
{
  ASSERT_LT_DBG( function_no, fe_values_.n_dofs() );
  ASSERT_LT_DBG( point_no, fe_values_.n_points() );
  arma::vec::fixed<spacedim> div;
  div.zeros();
  for (unsigned int c=0; c<spacedim*spacedim; ++c)
    div(c%spacedim) += fe_values_.shape_grad_component(function_no, point_no, first_tensor_component_+c)[first_tensor_component_+c/spacedim];
  return div;
}

template<unsigned int spacedim>
const FEValues<spacedim> &FEValuesViews::Tensor<spacedim>::base() const
{ return fe_values_; }
  




template class FEValuesViews::Scalar<3>;
template class FEValuesViews::Vector<3>;
template class FEValuesViews::Tensor<3>;


