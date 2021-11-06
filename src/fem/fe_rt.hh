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
 * @file    fe_rt.hh
 * @brief   Definitions of Raviart-Thomas finite elements.
 * @author  Jan Stebel
 */

#ifndef FE_RT_HH_
#define FE_RT_HH_

#include "fem/finite_element.hh"
#include "system/logger.hh"

/**
 * Space of Raviart-Thomas polynomials of order 0 (affine functions).
 * The basis functions are defined as
 * 
 *    x, x - e_1, ..., x - e_d,
 * 
 * where x is the space variable, e_i is the i-th canonical basis vector in R^d
 * and d is @p space_dim_.
 */
class RT0_space : public FunctionSpace {
public:
    
    RT0_space(unsigned int dim);
    
    double basis_value(unsigned int basis_index,
                       const arma::vec &point,
                       unsigned int comp_index
                       ) const override;
    
    const arma::vec basis_grad(unsigned int basis_index,
                               const arma::vec &point,
                               unsigned int comp_index
                              ) const override;

    unsigned int dim() const override { return this->space_dim_+1; }
};


/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim>
class FE_RT0 : public FiniteElement<dim>
{
public:

    /**
     * @brief Constructor.
     */
    FE_RT0();

};




/**
 * @brief Discontinuous Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim>
class FE_RT0_disc : public FiniteElement<dim>
{
public:

    /**
     * @brief Constructor.
     */
    FE_RT0_disc();

};



















#endif /* FE_RT_HH_ */
