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
 * @file    fe_p.hh
 * @brief   Definitions of basic Lagrangean finite elements with polynomial shape functions.
 * @author  Jan Stebel
 */

#ifndef FE_P_HH_
#define FE_P_HH_

#include <boost/exception/detail/error_info_impl.hpp>  // for error_info
#include <boost/exception/info.hpp>                    // for operator<<
#include <string>                                      // for string
#include <vector>                                      // for vector
#include <armadillo>
#include "fem/finite_element.hh"                       // for FiniteElement
#include "system/exc_common.hh"                        // for ExcAssertMsg
#include "system/exceptions.hh"                        // for ExcAssertMsg::...
#include "system/global_defs.h"                        // for OLD_ASSERT, msg

template <unsigned int dim, unsigned int spacedim> class FE_P0_XFEM;

/**
 * @brief Space of polynomial functions.
 *
 * This class serves for evaluation of the value and gradient
 * of a polynomial of order @p degree in @p dim variables.
 */
class PolynomialSpace : public FunctionSpace
{
public:

	/**
	 * @brief Constructor.
	 *
	 * Creates the coefficients of the basis.
	 */
    PolynomialSpace(unsigned int degree, unsigned int dim);

    const double basis_value(unsigned int basis_index,
                             const arma::vec &point,
                             unsigned int comp_index = 0
                            ) const override;
    
    const arma::vec basis_grad(unsigned int basis_index,
                               const arma::vec &point,
                               unsigned int comp_index = 0
                              ) const override;

    const unsigned int dim() const override { return powers.size(); }

private:

    /// Max. degree of polynomials.
    const unsigned int degree_;
    
    /**
     * @brief Coefficients of basis functions.
     *
     * Powers of x, y, z, ... in the i-th basis function are stored
     * in powers[i].
     */
    std::vector<arma::uvec> powers;

};






/**
 * @brief Conforming Lagrangean finite element on @p dim dimensional simplex.
 *
 * The finite element functions are continuous across the interfaces.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_P : public FiniteElement<dim,spacedim>
{
public:
    /// Constructor.
    FE_P(unsigned int degree);
    
protected:
    
    void init_dofs();

    /// Maximum degree of polynomials.
    unsigned int degree_;
    
    friend class FE_P0_XFEM<dim,spacedim>;
};


/**
 * @brief Discontinuous Lagrangean finite element on @p dim dimensional simplex.
 *
 * No continuity of the finite element functions across the interfaces is
 * imposed.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_P_disc : public FE_P<dim,spacedim>
{
public:

    /// Constructor.
    FE_P_disc(unsigned int degree);
    
    friend class FE_P0_XFEM<dim,spacedim>;
};














































#endif /* FE_P_HH_ */
