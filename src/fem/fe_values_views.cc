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

#include "fem/fe_values_views.hh"
#include "fem/fe_values.hh"

using namespace FEValuesViews;

template<unsigned int dim, unsigned int spacedim>
double FEValuesViews::Scalar<dim,spacedim>::value(unsigned int function_no, unsigned int point_no) const
{
  return fe_values_.shape_value_component(function_no, point_no, component_);
}

template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesViews::Scalar<dim,spacedim>::grad(unsigned int function_no, unsigned int point_no) const
{
  return fe_values_.shape_grad_component(function_no, point_no, component_);
}

template<unsigned int dim, unsigned int spacedim>
FEValuesBase<dim,spacedim> &FEValuesViews::Scalar<dim,spacedim>::base() const
{ return fe_values_; }
  



template<unsigned int dim, unsigned int spacedim>
arma::vec::fixed<spacedim> FEValuesViews::Vector<dim,spacedim>::value(unsigned int function_no, unsigned int point_no) const
{
  arma::vec::fixed<spacedim> v;
  for (unsigned int c=0; c<spacedim; ++c)
    v(c) = fe_values_.shape_value_component(function_no, point_no, first_vector_component_+c);
  return v;
}

template<unsigned int dim, unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Vector<dim,spacedim>::grad(unsigned int function_no, unsigned int point_no) const
{
  arma::mat::fixed<spacedim,spacedim> g;
  for (unsigned int c=0; c<spacedim; ++c)
    g.row(c) = fe_values_.shape_grad_component(function_no, point_no, first_vector_component_+c);
  return g;
}

template<unsigned int dim, unsigned int spacedim>
arma::mat::fixed<spacedim,spacedim> FEValuesViews::Vector<dim,spacedim>::sym_grad(unsigned int function_no, unsigned int point_no) const
{
  arma::mat::fixed<spacedim,spacedim> g = grad(function_no, point_no);
  return 0.5*(g+trans(g));
}

template<unsigned int dim, unsigned int spacedim>
double FEValuesViews::Vector<dim,spacedim>::divergence(unsigned int function_no, unsigned int point_no) const
{
  double div = 0;
  for (unsigned int c=0; c<spacedim; ++c)
    div += fe_values_.shape_grad_component(function_no, point_no, first_vector_component_+c)(first_vector_component_+c);
  return div;
}

template<unsigned int dim, unsigned int spacedim>
FEValuesBase<dim,spacedim> &FEValuesViews::Vector<dim,spacedim>::base() const
{ return fe_values_; }
  




template class FEValuesViews::Scalar<1,3>;
template class FEValuesViews::Scalar<2,3>;
template class FEValuesViews::Scalar<3,3>;

template class FEValuesViews::Vector<1,3>;
template class FEValuesViews::Vector<2,3>;
template class FEValuesViews::Vector<3,3>;


