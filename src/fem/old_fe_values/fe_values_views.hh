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
 * @file    fe_system.hh
 * @brief   Class FESystem for compound finite elements.
 * @author  Jan Stebel
 */

#ifndef FE_VALUES_VIEWS_HH
#define FE_VALUES_VIEWS_HH

#include <vector>
#include <armadillo>

template<unsigned int spacedim> class FEValues;



/**
 * FEValuesViews classes realize access to shape functions and their derivatives
 * for scalar, vector, tensor etc. part of the finite element.
 * 
 * The classes are motivated by deal.II.
 */
namespace FEValuesViews {

  template<unsigned int spacedim = 3>
  class Scalar {
    
  public:
    
    Scalar(const FEValues<spacedim> &fe_values, unsigned int component)
      : fe_values_(fe_values),
        component_(component)
    {};
    
    /**
     * @brief Return value of scalar shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    double value(unsigned int function_no, unsigned int point_no) const;
    
    /**
     * @brief Return gradient of scalar shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::vec::fixed<spacedim> grad(unsigned int function_no, unsigned int point_no) const;
    
    /// Returns the FEValues class.
    const FEValues<spacedim> &base() const;
    
  private:
    
    /// Base FEValues class for access to the FE.
    const FEValues<spacedim> &fe_values_;
    
    /// Index of the scalar component.
    unsigned int component_;
  };


  template<unsigned int spacedim = 3>
  class Vector {
    
  public:
    
    Vector(const FEValues<spacedim> &fe_values, unsigned int component)
      : fe_values_(fe_values),
        first_vector_component_(component)
    {};
    
    /**
     * @brief Return value of vector-valued shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::vec::fixed<spacedim> value(unsigned int function_no, unsigned int point_no) const;
    
    /**
     * @brief Return gradient of vector-valued shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::mat::fixed<spacedim,spacedim> grad(unsigned int function_no, unsigned int point_no) const;
    
    /**
     * @brief Return symmetric gradient of vector-valued shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::mat::fixed<spacedim,spacedim> sym_grad(unsigned int function_no, unsigned int point_no) const;

    /**
     * @brief Return divergence of vector-valued shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    double divergence(unsigned int function_no, unsigned int point_no) const;
    
    /// Returns the FEValues class.
    const FEValues<spacedim> &base() const;
    
  private:
    
    /// Base FEValues class for access to the FE.
    const FEValues<spacedim> &fe_values_;
    
    /// Index of the first component of the vector.
    unsigned int first_vector_component_;
  };
  
  
  template<unsigned int spacedim = 3>
  class Tensor {
      
  public:
    
    Tensor(const FEValues<spacedim> &fe_values, unsigned int component)
      : fe_values_(fe_values),
        first_tensor_component_(component)
    {};
    
    /**
     * @brief Return value of tensor-valued shape function.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::mat::fixed<spacedim,spacedim> value(unsigned int function_no, unsigned int point_no) const;
    
    /**
     * @brief Return partial derivative of tensor-valued shape function.
     * @param variable_no Index of spacial variable w.r. to which we differentiate.
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::mat::fixed<spacedim,spacedim> derivative(
        unsigned int variable_no,
        unsigned int function_no,
        unsigned int point_no) const;
    
    /**
     * @brief Return divergence of tensor-valued shape function.
     * 
     * The result is a vector whose components are divergences of tensor columns, i.e.
     * (div T)_i = dT_ji / dx_j.
     * 
     * @param function_no Index of shape function within the FE.
     * @param point_no    Index of quadrature point.
     */
    arma::vec::fixed<spacedim> divergence(unsigned int function_no, unsigned int point_no) const;
    
    /// Returns the FEValues class.
    const FEValues<spacedim> &base() const;
    
  private:
    
    /// Base FEValues class for access to the FE.
    const FEValues<spacedim> &fe_values_;
    
    /// Index of the first component of the vector.
    unsigned int first_tensor_component_;
    
  };

}






#endif
