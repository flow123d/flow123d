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


class RT0_space : public FunctionSpace {
public:
    
    RT0_space(unsigned int dim);
    
    const double basis_value(unsigned int basis_index,
                             const arma::vec &point,
                             unsigned int comp_index
                            ) const override;
    
    const arma::vec basis_grad(unsigned int basis_index,
                               const arma::vec &point,
                               unsigned int comp_index
                              ) const override;

    const unsigned int dim() const override { return this->space_dim_+1; }
};


/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0 : public FiniteElement<dim,spacedim>
{
public:

    /// Number of raw basis functions.
    static const unsigned int n_raw_functions = dim+1;

    /**
     * @brief Constructor.
     */
    FE_RT0();

    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Destructor.
     */
    ~FE_RT0() override {}
    
private:
  
    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param q Quadrature.
     * @param data The precomputed finite element data on the reference cell.
     * @param fv_data The data to be computed.
     */
    void fill_fe_values(const Quadrature<dim> &q,
            FEInternalData &data,
            FEValuesData<dim,spacedim> &fv_data) override;

    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &q) override;
};
























#endif /* FE_RT_HH_ */
