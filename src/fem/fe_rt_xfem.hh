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
 * @file    fe_rt_xfem.hh
 * @brief   Definitions of enriched Raviart-Thomas finite elements.
 * @author  Pavel Exner
 */

#ifndef FE_RT_XFEM_HH_
#define FE_RT_XFEM_HH_

#include "fem/finite_element.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_p.hh"

#include "fem/global_enrichment_func.hh"
#include "fem/fe_values.hh"
#include "quadrature/qxfem.hh"
#include "mesh/ref_element.hh"

#include "system/logger.hh"

template <unsigned int dim, unsigned int spacedim> class FE_RT0_XFEM;

/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0_XFEM : public FiniteElement<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
//     using FiniteElement<dim,spacedim>::number_of_pairs;
//     using FiniteElement<dim,spacedim>::number_of_triples;
//     using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::generalized_support_points;
    using FiniteElement<dim,spacedim>::order;
    using FiniteElement<dim,spacedim>::is_scalar_fe;
//     using FiniteElement<dim,spacedim>::node_matrix;
    
    
    FE_RT0<dim,spacedim> rt0;
    FE_P_disc<1,dim, spacedim> pu;
    
    std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr;
    
    unsigned int n_regular_dofs;
    
public:

    /// Number of raw basis functions.
    static const unsigned int n_raw_functions = dim+1;

    /**
     * @brief Constructor.
     */
    FE_RT0_XFEM(std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr);

    /**
     * @brief The scalar variant of basis_vector must be implemented but may not be used.
     */
    double basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief The scalar variant of basis_grad_vector must be implemented but may not be used.
     */
    arma::vec::fixed<dim> basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the @p ith basis function evaluated at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::vec::fixed<dim> basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the gradient of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const;

    /**
     * @brief Returns the divergence of the @p ith basis function at the point @p p.
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    double basis_div(const unsigned int i, const arma::vec::fixed<dim> &p) const;
    
    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature.
     * @param flags Flags that indicate what quantities should be calculated.
     */
    FEInternalData *initialize(const Quadrature<dim> &quad, UpdateFlags flags);

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
    virtual void fill_fe_values(const QXFEM<dim,spacedim> &quad,
            FEInternalData &data,
            FEValuesData<dim,spacedim> &fv_data);

    /**
     * @brief Destructor.
     */
    virtual ~FE_RT0_XFEM() {};

};










template <unsigned int dim, unsigned int spacedim>
FE_RT0_XFEM<dim,spacedim>::FE_RT0_XFEM(std::vector<GlobalEnrichmentFunc<dim,spacedim>*> enr)
: enr(enr)
{
//     arma::vec::fixed<dim> sp;

    this->init();

    n_regular_dofs = dim+1;
    
    // regular + enriched from every singularity
    number_of_dofs = n_regular_dofs + n_regular_dofs * enr.size();
    number_of_single_dofs[dim] = number_of_dofs;

//     for (unsigned int sid=0; sid<RefElement<dim>::n_sides; ++sid)
//     {
//         sp.fill(0);
//         for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; ++i)
//             sp += RefElement<dim>::node_coords(RefElement<dim>::template interact<0,dim-1>(sid)[i]);
//         sp /= RefElement<dim>::n_nodes_per_side;
//         generalized_support_points.push_back(sp);
//     }

    generalized_support_points = rt0.generalized_support_points;
    order = 1;
    is_scalar_fe = false;
}

template <unsigned int dim, unsigned int spacedim>
double FE_RT0_XFEM<dim,spacedim>::basis_value(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT_DBG(false).error("basis_value() may not be called for vectorial finite element.");

    return 0.0;
}

template <unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0_XFEM<dim,spacedim>::basis_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT_DBG(false).error("basis_grad() may not be called for vectorial finite element.");
    return arma::vec::fixed<dim>();
}

template <unsigned int dim, unsigned int spacedim>
arma::vec::fixed<dim> FE_RT0_XFEM<dim,spacedim>::basis_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT_DBG(i<n_raw_functions).error("Index of basis function is out of range.");

    arma::vec::fixed<dim> v(p);
    
    if (i > 0)
        v[i-1] -= 1;

    return v;
}

template <unsigned int dim, unsigned int spacedim>
arma::mat::fixed<dim,dim> FE_RT0_XFEM<dim,spacedim>::basis_grad_vector(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT_DBG(i<n_raw_functions).error("Index of basis function is out of range.");

    return arma::eye(dim,dim);
}

template <unsigned int dim, unsigned int spacedim>
double FE_RT0_XFEM<dim,spacedim>::basis_div(const unsigned int i, const arma::vec::fixed<dim> &p) const
{
    ASSERT_DBG(i<n_raw_functions).error("Index of basis function is out of range.");

    return dim;
}


template <unsigned int dim, unsigned int spacedim>
FEInternalData *FE_RT0_XFEM<dim,spacedim>::initialize(const Quadrature<dim> &quad, UpdateFlags flags)
{
    ASSERT_DBG(false).error("No internal data on reference element for XFEM.");
//     FEInternalData *data = new FEInternalData;
// 
//     //TODO: fill values directly without these objects..
//     FEInternalData *data_rt = rt0.initialize(quad,flags);
//     FEInternalData *data_pu = pu.initialize(quad,flags);
//     
//     data->basis_vectors = data_rt->basis_vectors;
//     data->basis_grad_vectors = data_rt->basis_grad_vectors;
//     data->basis_values = data_pu->basis_values;
//     data->basis_grads = data_pu->basis_grads;
//     
//     delete data_rt;
//     delete data_pu;
//     
//     return data;
    return nullptr;
}

template <unsigned int dim, unsigned int spacedim>
inline UpdateFlags FE_RT0_XFEM<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_values)
        f |= update_jacobians | update_volume_elements;

    if (flags & update_gradients)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements;

    return f;
}

template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM<dim,spacedim>::fill_fe_values(
        const QXFEM<dim,spacedim> &quad,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    ElementFullIter ele = *fv_data.present_cell;
    typedef typename Space<spacedim>::Point Point;
    unsigned int j;
    
    // can we suppose for this FE and element that:
    //  - jacobian (and its inverse and determinant) is constant on the element
    //  - abuse mapping to compute the normals
    
    
    // shape values
    if (fv_data.update_flags & update_values)
    {
        // compute normals - TODO: this should do mapping, but it does for fe_side_values..
        vector<arma::vec::fixed<spacedim>> normals(RefElement<dim>::n_sides);
        for (j = 0; j < RefElement<dim>::n_sides; j++){
            normals[j] = trans(fv_data.inverse_jacobians[0])*RefElement<dim>::normal_vector(j);
            normals[j] = normals[j]/norm(normals[j],2);
        }
        
        // for SGFEM
        // values of enrichment function at generalized_support_points
        // here: n_regular_dofs = RefElement<dim>::n_nodes = RefElement<dim>::n_sides = generalized_support_points.size()
        vector<vector<double>> enr_dof_val(enr.size());
        for (unsigned int w=0; w<enr.size(); w++){
            enr_dof_val.resize(RefElement<dim>::n_sides);
            for (unsigned int i = 0; i < RefElement<dim>::n_sides; i++)
            {
                Point real_point;
                real_point.zeros();
                for (j = 0; j < n_regular_dofs; j++)
                    real_point = real_point + rt0.generalized_support_points[j] * ele->node[j]->point();
                
                enr_dof_val[w][i] = arma::dot(enr[w]->vector(real_point), normals[i]);
            }
        }

        
        vector<arma::vec::fixed<spacedim> > vectors(number_of_dofs);
        
        for (unsigned int q = 0; q < quad.size(); q++)
        {
            // compute PU
            arma::vec pu_values(RefElement<dim>::n_nodes);
            for (j=0; j<RefElement<dim>::n_nodes; j++)
                    pu_values[j] = pu.basis_value(j, quad.point(q));
            pu_values = pu.node_matrix * pu_values;
            
            // compute PU grads
            arma::mat pu_grads(RefElement<dim>::n_nodes, dim);
            for (j=0; j<RefElement<dim>::n_nodes; j++)
                pu_grads.row(j) = arma::trans(pu.basis_grad(j, quad.point(q)));
            pu_grads = pu.node_matrix * pu_grads;
            // real_pu_grad = pu_grads[i] * fv_data.inverse_jacobians[i];
            
            
            //fill regular shape functions
            for (j=0; j<n_regular_dofs; j++)
                vectors[j] = fv_data.jacobians[q] * rt0.basis_vector(j,quad.point(q)) / fv_data.determinants[q];
            
            j = n_regular_dofs;
            for (unsigned int w=0; w<enr.size(); w++)
            {
                //compute interpolant
                arma::vec::fixed<spacedim> interpolant;
                interpolant.zeros();
                for (unsigned int k=0; k<n_regular_dofs; k++)
                    interpolant += vectors[k] * enr_dof_val[w][k];
                
                for (unsigned int k=0; k<n_regular_dofs; k++)
                {
                    j++;
                    vectors[j] =  pu_values[k] * (enr[w]->value(quad.real_point(q)) - interpolant);
                }
            }   
                
            fv_data.shape_vectors[q] = vectors;
        }
    }

//     // shape gradients
//     if (fv_data.update_flags & update_gradients)
//     {
//         vector<arma::mat::fixed<spacedim,spacedim> > grads;
//         grads.resize(dim+1);
//         for (unsigned int q = 0; q < quad.size(); q++)
//         {
//             for (unsigned int k=0; k<dim+1; k++)
//                 grads[k] = fv_data.jacobians[q]*data.basis_grad_vectors[q][k]*fv_data.inverse_jacobians[q]/fv_data.determinants[q];
// 
//             fv_data.shape_grad_vectors[q] = grads;
//         }
//     }
}













#endif // FE_RT_XFEM_HH_
