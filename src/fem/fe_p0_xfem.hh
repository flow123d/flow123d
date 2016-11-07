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

#ifndef FE_P0_XFEM_HH_
#define FE_P0_XFEM_HH_

#include "fem/finite_element_enriched.hh"
#include "fem/fe_p.hh"

#include "fem/global_enrichment_func.hh"
#include "fem/fe_values.hh"
#include "quadrature/qxfem.hh"
#include "mesh/ref_element.hh"

#include "system/logger.hh"

template <unsigned int dim, unsigned int spacedim> class FE_P0_XFEM;

/**
 * @brief P0 element of order 0.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_P0_XFEM : public FiniteElementEnriched<dim,spacedim>
{
    using FiniteElement<dim,spacedim>::number_of_dofs;
    using FiniteElement<dim,spacedim>::number_of_single_dofs;
//     using FiniteElement<dim,spacedim>::number_of_pairs;
//     using FiniteElement<dim,spacedim>::number_of_triples;
//     using FiniteElement<dim,spacedim>::number_of_sextuples;
    using FiniteElement<dim,spacedim>::generalized_support_points;
    using FiniteElement<dim,spacedim>::unit_support_points;
    using FiniteElement<dim,spacedim>::order;
    using FiniteElement<dim,spacedim>::is_scalar_fe;
//     using FiniteElement<dim,spacedim>::node_matrix;
    
    using FiniteElementEnriched<dim,spacedim>::fe;
    using FiniteElementEnriched<dim,spacedim>::pu;
    using FiniteElementEnriched<dim,spacedim>::enr;
    using FiniteElementEnriched<dim,spacedim>::n_regular_dofs_;
    
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
public:

    /**
     * @brief Constructor.
     */
    FE_P0_XFEM(FE_P_disc<0,dim,spacedim>* fe,std::vector<EnrichmentPtr> enr);


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
    void fill_fe_values(const Quadrature<dim> &quad,
                        FEInternalData &data,
                        FEValuesData<dim,spacedim> &fv_data) override;

    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param q Quadrature.
     * @param data The precomputed finite element data on the reference cell.
     * @param fv_data The data to be computed.
     */
    void fill_fe_values(const Quadrature<dim> &quad,
                        FEValuesData<dim,spacedim> &fv_data);
    
    /**
     * @brief Destructor.
     */
    ~FE_P0_XFEM() {};

};






template <unsigned int dim, unsigned int spacedim>
FE_P0_XFEM<dim,spacedim>::FE_P0_XFEM(FE_P_disc<0,dim,spacedim>* fe,
                                     std::vector<EnrichmentPtr> enr)
: FiniteElementEnriched<dim,spacedim>(fe,enr)
{
    order = 0;
    is_scalar_fe = true;
}


template <unsigned int dim, unsigned int spacedim>
inline UpdateFlags FE_P0_XFEM<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_values)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    if (flags & update_gradients)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    return f;
}


template <unsigned int dim, unsigned int spacedim>
inline void FE_P0_XFEM<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    fill_fe_values(quad,fv_data);
}

template <unsigned int dim, unsigned int spacedim>
inline void FE_P0_XFEM<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEValuesData<dim,spacedim> &fv_data)
{
    ElementFullIter ele = *fv_data.present_cell;
    typedef typename Space<spacedim>::Point Point;
    unsigned int j;
    
    // can we suppose for this FE and element that:
    //  - jacobian (and its inverse and determinant) is constant on the element
    //  - abuse mapping to compute the normals
    
//     DBGCOUT("FE_P0_XFEM fill fe_values\n");
    
    // shape values
    if (fv_data.update_flags & update_values)
    {
//         DBGMSG("interpolation\n");
        // for SGFEM
        // values of enrichment function at generalized_support_points
        // here: n_regular_dofs = RefElement<dim>::n_nodes = RefElement<dim>::n_sides = generalized_support_points.size()
        vector<vector<double>> enr_dof_val(enr.size());
        auto& gen_points = fe->get_generalized_support_points();
        for (unsigned int w=0; w<enr.size(); w++){
            enr_dof_val[w].resize(gen_points.size()); // for P0 is equal 1 (barycenter)
            for (unsigned int i = 0; i < gen_points.size(); i++)
            {
                // compute real generalized_support_point
                Point real_point; real_point.zeros();
                
                arma::vec::fixed<dim+1> bp = RefElement<dim>::local_to_bary(gen_points[i]);
                for (j = 0; j < RefElement<dim>::n_nodes; j++)
                    real_point += bp[j] * ele->node[j]->point();
                    
                enr_dof_val[w][i] = enr[w]->value(real_point);
//                 cout << "interpolant: [" << i << "]: " << setprecision(15) << enr_dof_val[w][i] << endl;
            }
        }
        
        vector<double> values(number_of_dofs);
        
        for (unsigned int q = 0; q < quad.size(); q++)
        {
//             DBGMSG("pu q[%d]\n",q);
            // compute PU
            arma::vec pu_values(pu.n_dofs());
            for (j=0; j < pu.n_dofs(); j++)
                pu_values[j] = pu.basis_value(j, quad.point(q));
            pu_values = pu.get_node_matrix() * pu_values;
            
//             DBGMSG("pu grad q[%d]\n",q);
            // compute PU grads
//             arma::mat pu_grads(RefElement<dim>::n_nodes, dim);
//             for (j=0; j<RefElement<dim>::n_nodes; j++)
//                 pu_grads.row(j) = arma::trans(pu.basis_grad(j, quad.point(q)));
//             pu_grads = pu.node_matrix * pu_grads;
            // real_pu_grad = pu_grads[i] * fv_data.inverse_jacobians[i];
            
            
//             DBGMSG("regular shape vectors q[%d]\n",q);
            //fill regular shape functions
            for (j=0; j<n_regular_dofs_; j++)
                values[j] = fe->basis_value(j,quad.point(q));
            
            j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
//                 DBGMSG("interpolant w[%d] q[%d]\n",w,q);
                //compute interpolant
                double interpolant = 0.0;
                for (unsigned int k=0; k < n_regular_dofs_; k++)
                    interpolant += values[k] * enr_dof_val[w][k];
                
//                 quad.real_point(q).print(cout);
//                 DBGMSG("enriched shape value w[%d] q[%d]\n",w,q);
                for (unsigned int k=0; k < pu.n_dofs(); k++)
                {
                    values[j] +=  pu_values[k] * (enr[w]->value(fv_data.points[q]) - interpolant);
//                     values[j] =  pu_values[k];
//                     values[j] =  interpolant;//(enr[w]->value(fv_data.points[q]) - interpolant);
//                     values[j] =  enr[w]->value(fv_data.points[q]);
                    j++;
                }
            }   
                
            fv_data.shape_values[q] = values;
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

#endif // FE_P0_XFEM_HH_
