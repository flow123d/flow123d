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
 * @file    finite_element_enriched.hh
 * @brief   Definitions of common enriched finite element.
 * @author  Pavel Exner
 */

#ifndef FINITE_ELEMENT_ENRICHED_HH_
#define FINITE_ELEMENT_ENRICHED_HH_

#include "fem/finite_element.hh"
#include "fem/global_enrichment_func.hh"
#include "fem/fe_p.hh"

template <unsigned int dim, unsigned int spacedim> class FiniteElementEnriched;


template <unsigned int dim, unsigned int spacedim>
class FiniteElementEnriched : public FiniteElement<dim,spacedim>
{
protected:
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
    
    FiniteElement<dim,spacedim> *fe;
    FE_P_disc<1,dim, spacedim> pu;
    
    unsigned int n_regular_dofs_;
    
    std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr;
public:
    /**
     * @brief Constructor.
     */
    FiniteElementEnriched(FiniteElement<dim,spacedim>* fe,
                          std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr);
    
    virtual ~FiniteElementEnriched(){}
    
    const unsigned int n_regular_dofs() const;
    const unsigned int n_enriched_dofs() const;
    
    
    /**
     * @brief The scalar variant of basis_vector must be implemented but may not be used.
     */
    double basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const override;

    /**
     * @brief The scalar variant of basis_grad_vector must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const override;

    /**
     * @brief Returns the @p ith basis function evaluated at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::vec::fixed<dim> basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const override;

    /**
     * @brief Returns the gradient of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const override;

    /**
     * @brief Returns the divergence of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    double basis_div(const unsigned int i, const arma::vec::fixed<dim> &p) const override;
    
    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &quad, UpdateFlags flags) override;
};

template <unsigned int dim, unsigned int spacedim>
FiniteElementEnriched<dim,spacedim>::FiniteElementEnriched(FiniteElement<dim,spacedim>* fe,
                                                           std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr)
: fe(fe), enr(enr)
{
    this->init();

    n_regular_dofs_ = fe->n_dofs();
    
    // regular + enriched from every singularity
    number_of_dofs = n_regular_dofs_ + pu.n_dofs() * enr.size();
    number_of_single_dofs[dim] = number_of_dofs;
}

template <unsigned int dim, unsigned int spacedim>
const unsigned int FiniteElementEnriched<dim,spacedim>::n_regular_dofs() const
{   return n_regular_dofs_;}

template <unsigned int dim, unsigned int spacedim>
const unsigned int FiniteElementEnriched<dim,spacedim>::n_enriched_dofs() const
{   return number_of_dofs - n_regular_dofs_;}

template <unsigned int dim, unsigned int spacedim>
double FiniteElementEnriched<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    return fe->basis_value(i,p);
}

template <unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FiniteElementEnriched<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    return fe->basis_grad(i,p);
}

template <unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FiniteElementEnriched<dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    return fe->basis_vector(i,p);
}

template <unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FiniteElementEnriched<dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    return fe->basis_grad_vector(i,p);
}

template <unsigned int dim, unsigned int spacedim>
double FiniteElementEnriched<dim,spacedim>::basis_div(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    return fe->basis_div(i,p);
}

template <unsigned int dim, unsigned int spacedim>
FEInternalData *FiniteElementEnriched<dim,spacedim>::initialize(const Quadrature<dim> &quad, UpdateFlags flags)
{
//     ASSERT_DBG(false).error("No internal data on reference element for XFEM.");
    // return empty internal data since all internal data are element dependent
    return new FEInternalData;
}


#endif // FINITE_ELEMENT_ENRICHED_HH_