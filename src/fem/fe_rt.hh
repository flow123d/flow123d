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
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0 : public FiniteElement<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
    using FiniteElement<dim,spacedim>::number_of_pairs;
    using FiniteElement<dim,spacedim>::number_of_triples;
    using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::generalized_support_points;
    using FiniteElement<dim,spacedim>::is_primitive_;
    using FiniteElement<dim,spacedim>::node_matrix;

public:

    /// Number of raw basis functions.
    static const unsigned int n_raw_functions = dim+1;

    /**
     * @brief Constructor.
     */
    FE_RT0();

    /**
     * @brief Calculates the value of the @p comp-th component of
     * the @p i-th raw basis function at the
     * point @p p on the reference element (for vector-valued FE).
     *
     * @param i    Number of the basis function.
     * @param p    Point of evaluation.
     * @param comp Number of vector component.
     */
    double basis_value(const unsigned int i,
            const arma::vec::fixed<dim> &p, const unsigned int comp) const override;

    /**
     * @brief Calculates the @p comp-th component of the gradient
     * of the @p i-th raw basis function at the point @p p on the
     * reference element (for vector-valued FE).
     *
     * @param i    Number of the basis function.
     * @param p    Point of evaluation.
     * @param comp Number of vector component.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i,
            const arma::vec::fixed<dim> &p, const unsigned int comp) const override;

    /**
     * @brief Computes the conversion matrix from internal basis to shape functions.
     *
     * Initializes the @p node_matrix for computing the coefficients
     * of the raw basis functions from values at support points.
     * Since Raviart-Thomas element is not Lagrangean, the method has
     * to be reimplemented.
     */
    void compute_node_matrix();

    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &q) override;
    
    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param q Quadrature.
     * @param data The precomputed finite element data on the reference cell.
     * @param fv_data The data to be computed.
     */
    virtual void fill_fe_values(const Quadrature<dim> &q,
            FEInternalData &data,
            FEValuesData<dim,spacedim> &fv_data);

    /**
     * @brief Destructor.
     */
    virtual ~FE_RT0() {};

};
























#endif /* FE_RT_HH_ */
