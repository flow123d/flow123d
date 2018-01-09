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
    
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
    FiniteElement<dim,spacedim> *fe;
//     FE_P_disc<1,dim, spacedim> pu;
    FE_P_disc<0,dim, spacedim> pu;
    
    unsigned int n_regular_dofs_;
    
    std::vector<EnrichmentPtr> enr;
public:
    /**
     * @brief Constructor.
     */
    FiniteElementEnriched(FiniteElement<dim,spacedim>* fe,
                          std::vector<EnrichmentPtr> enr);
    
    virtual ~FiniteElementEnriched(){}
    
    const unsigned int n_regular_dofs() const;
    const unsigned int n_enriched_dofs() const;
    
    
    /**
     * @brief The scalar variant of basis_vector must be implemented but may not be used.
     */
    double basis_value(const unsigned int i,
                       const arma::vec::fixed<dim> &p,
                       const unsigned int comp = 0) const override;

    /**
     * @brief The scalar variant of basis_grad_vector must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i,
                                     const arma::vec::fixed<dim> &p,
                                     const unsigned int comp = 0) const override;

    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &quad) override;
};

template <unsigned int dim, unsigned int spacedim>
FiniteElementEnriched<dim,spacedim>::FiniteElementEnriched(FiniteElement<dim,spacedim>* fe,
                                                           std::vector<EnrichmentPtr> enr)
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
double FiniteElementEnriched<dim,spacedim>::basis_value(const unsigned int i,
                                                        const arma::vec::fixed<dim> &p,
                                                        const unsigned int comp) const
{
    return fe->basis_value(i,p,comp);
}

template <unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FiniteElementEnriched<dim,spacedim>::basis_grad(const unsigned int i,
                                                                      const arma::vec::fixed<dim> &p,
                                                                      const unsigned int comp) const
{
    return fe->basis_grad(i,p,comp);
}

template <unsigned int dim, unsigned int spacedim>
FEInternalData *FiniteElementEnriched<dim,spacedim>::initialize(const Quadrature<dim> &quad)
{
//     ASSERT_DBG(false).error("No internal data on reference element for XFEM.");
    // return empty internal data since all internal data are element dependent
    return new FEInternalData;
}


#endif // FINITE_ELEMENT_ENRICHED_HH_
