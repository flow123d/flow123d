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

#ifndef FE_RT0_XFEM_HH_
#define FE_RT0_XFEM_HH_

#include "fem/finite_element_enriched.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_p.hh"

#include "fem/global_enrichment_func.hh"
#include "fem/fe_values.hh"
#include "quadrature/qxfem.hh"
#include "mesh/ref_element.hh"

#include "system/logger.hh"

#include "quadrature/quadrature_lib.hh"
#include "fem/mapping_p1.hh"

template <unsigned int dim, unsigned int spacedim> class FE_RT0_XFEM;

/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0_XFEM : public FiniteElementEnriched<dim,spacedim>
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
    
    using FiniteElementEnriched<dim,spacedim>::fe;
    using FiniteElementEnriched<dim,spacedim>::pu;
    using FiniteElementEnriched<dim,spacedim>::enr;
    using FiniteElementEnriched<dim,spacedim>::n_regular_dofs_;
    
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
public:

    /**
     * @brief Constructor.
     */
    FE_RT0_XFEM(FE_RT0<dim,spacedim>* fe,std::vector<EnrichmentPtr> enr);


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
    ~FE_RT0_XFEM() {};

};






template <unsigned int dim, unsigned int spacedim>
FE_RT0_XFEM<dim,spacedim>::FE_RT0_XFEM(FE_RT0<dim,spacedim>* fe,
                                       std::vector<EnrichmentPtr> enr)
: FiniteElementEnriched<dim,spacedim>(fe,enr)
{
    order = 1;
    is_scalar_fe = false;
}


template <unsigned int dim, unsigned int spacedim>
inline UpdateFlags FE_RT0_XFEM<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_values)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    if (flags & update_gradients)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    return f;
}


template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    fill_fe_values(quad,fv_data);
}

template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEValuesData<dim,spacedim> &fv_data)
{
    ElementFullIter ele = *fv_data.present_cell;
    typedef typename Space<spacedim>::Point Point;
    
    // can we suppose for this FE and element that:
    //  - jacobian (and its inverse and determinant) is constant on the element
    //  - abuse mapping to compute the normals
    
//     DBGCOUT("FE_RT0_XFEM fill fe_values\n");
    
//     pu.get_node_matrix().print(cout,"pu_node_matrix");
//     fe->get_node_matrix().print(cout,"fe_rt_node_matrix");
    
    
    vector<arma::vec::fixed<spacedim>> normals;
    vector<vector<double>> enr_dof_val;
    
    if (fv_data.update_flags & (update_values | update_divergence | update_gradients))
    {
        normals.resize(RefElement<dim>::n_sides);
//         DBGCOUT("normals\n");
        // compute normals - TODO: this should do mapping, but it does it for fe_side_values..
        for (unsigned int j = 0; j < RefElement<dim>::n_sides; j++){
            normals[j] = trans(fv_data.inverse_jacobians[0])*RefElement<dim>::normal_vector(j);
            normals[j] = normals[j]/norm(normals[j],2);
//             normals[j].print(cout, "internal normal");
        }
        
        enr_dof_val.resize(enr.size());
        // for SGFEM we need to compute fluxes of glob. enr. function over faces
        ASSERT(dim == 2);
        //TODO :
        //  - create quadrature factory for uniform quad in 1d,2d,3d
        //  - create quadrature<dim> from quadrature<dim-1> (see mapping.transform_subquadrature() )
        QGauss<1> qside(1);
        const uint qsize=100;
        qside.resize(qsize);
        double qweight = 1.0/qsize;
        for(unsigned int q=0; q<qsize; q++){
            qside.set_point(q, arma::vec({0.5*qweight + q*qweight}));
            qside.set_weight(q,qweight);
        }
        MappingP1<2,3> map;
        
        for (unsigned int w=0; w<enr.size(); w++){
            enr_dof_val[w].resize(RefElement<dim>::n_sides);
            for (unsigned int j=0; j<RefElement<dim>::n_sides; j++){
                
                QGauss<2> quad_side(1);
                map.transform_subquadrature(j,ele->permutation_idx_[j],qside,quad_side);
                arma::mat ele_mat = map.element_map(*ele);
//                 ele_mat.print(cout,"ele_map");
//                 for(unsigned int q=0; q < quad_side.size(); q++){
//                     arma::vec qp = quad_side.point(q);
//                     cout << qp(0) << " " << qp(1) << " " << "\n";
//                 }
                    
                double side_measure = ele->side(j)->measure();
                double val = 0;
//                 DBGCOUT("edge quad [" << j << "]\n");
                for(unsigned int q=0; q < quad_side.size(); q++){
                    
                    arma::vec3 real_point = ele_mat.cols(1,dim) * quad_side.point(q) + ele_mat.col(0);
//                     cout << real_point(0) << " " << real_point(1) << " " << real_point(2) << "\n";
                    val += arma::dot(enr[w]->vector(real_point),normals[j])
                        // this makes JxW on the triangle side:
                        * quad_side.weight(q)
                        * side_measure;
                }
                enr_dof_val[w][j] = val;
//                 DBGVAR(val);
            }
        }
      
    }
//         //         DBGMSG("interpolation\n");
//         // for SGFEM
//         // values of enrichment function at generalized_support_points
//         // here: n_regular_dofs = RefElement<dim>::n_nodes = RefElement<dim>::n_sides = generalized_support_points.size()
//         auto& gen_points = fe->get_generalized_support_points();
//         for (unsigned int w=0; w<enr.size(); w++){
//             enr_dof_val[w].resize(gen_points.size());
//             for (unsigned int i = 0; i < gen_points.size(); i++)
//             {
//                 // compute real generalized_support_point
//                 Point real_point; real_point.zeros();
//                 
// //                 arma::vec::fixed<dim+1> bp = RefElement<dim>::local_to_bary(gen_points[i]);
// //                 //bp.print(cout,"bp");
// //                 for (j = 0; j < RefElement<dim>::n_nodes; j++)
// //                     real_point += bp[(j+1)%dim] * ele->node[j]->point();
//                 
//                 // compute barycenter of the side ( = real generalized_support_point)
//                 for (j = 0; j < RefElement<dim>::n_nodes_per_side; j++)
//                     real_point = real_point + ele->node[RefElement<dim>::interact(Interaction<0,1>(i))[j]]->point();
//                 real_point = real_point / RefElement<dim>::n_nodes_per_side;
//                 
// //                 real_point.print(cout,"real");
// //                 enr[w]->vector(real_point).print(cout);
// //                 normals[i].print(cout);
//                 enr_dof_val[w][i] = arma::dot(enr[w]->vector(real_point), normals[i]);
// //                 cout << "interpolant: [" << i << "]: " << setprecision(15) << enr_dof_val[w][i] << endl;
//             }
//         }
//     }
    
    
    // shape values
    if (fv_data.update_flags & update_values)
    {
        vector<arma::vec::fixed<spacedim> > vectors(number_of_dofs);
        
        for (unsigned int q = 0; q < quad.size(); q++)
        {
//             DBGMSG("pu q[%d]\n",q);
            // compute PU: PU=1
            
            
//             DBGMSG("regular shape vectors q[%d]\n",q);
            //fill regular shape functions
            arma::mat::fixed<dim+1,dim> raw_values;
            arma::mat::fixed<dim+1,dim> shape_values;
            
            for (unsigned int j=0; j<n_regular_dofs_; j++)
                raw_values.row(j) = trans(fe->basis_vector(j, quad.point(q)));
            
            shape_values = fe->get_node_matrix() * raw_values;
            
            unsigned int j;
            for (j=0; j<n_regular_dofs_; j++)
                vectors[j] = fv_data.jacobians[q] * trans(shape_values.row(j)) / fv_data.determinants[q];
            
            //fill enriched shape functions
            j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
//                 DBGMSG("interpolant w[%d] q[%d]\n",w,q);
                //compute interpolant
                arma::vec::fixed<spacedim> interpolant;
                interpolant.zeros();
                for (unsigned int k=0; k < n_regular_dofs_; k++)
                    interpolant += vectors[k] * enr_dof_val[w][k];
                
//                 quad.real_point(q).print(cout);
//                 DBGMSG("enriched shape value w[%d] q[%d]\n",w,q);
                for (unsigned int k=0; k < pu.n_dofs(); k++)
                {
                    vectors[j] =  enr[w]->vector(fv_data.points[q]) - interpolant;
//                     vectors[j] =  pu_values[k] * (enr[w]->vector(fv_data.points[q]) - interpolant);
//                     vectors[j] =  pu_values[k] * (enr[w]->vector(fv_data.points[q]));
//                     vectors[j] =  interpolant;//(enr[w]->vector(fv_data.points[q]) - interpolant);
//                     vectors[j] =  enr[w]->vector(fv_data.points[q]);
                    j++;
                }
            }   
                
            fv_data.shape_vectors[q] = vectors;
        }
    }

    // divergence
    if (fv_data.update_flags & update_divergence)
    {
//         arma::mat::fixed<dim,dim> unit_grad;
        vector<double> divs(number_of_dofs);
        
//         arma::vec pu_values(RefElement<dim>::n_nodes);
//         arma::mat pu_grads(RefElement<dim>::n_nodes, dim);
//         arma::mat real_pu_grads(RefElement<dim>::n_nodes, spacedim);
            
        for (unsigned int q = 0; q < quad.size(); q++)
        {
//             DBGMSG("pu q[%d]\n",q);
            // compute PU
//             for (unsigned int j=0; j<RefElement<dim>::n_nodes; j++)
//                     pu_values[j] = pu.basis_value(j, quad.point(q));
//             pu_values = pu.get_node_matrix() * pu_values;
//             
// //             DBGMSG("pu grad q[%d]\n",q);
//             // compute PU grads
//             for (unsigned int j=0; j<RefElement<dim>::n_nodes; j++)
//                 pu_grads.row(j) = arma::trans(pu.basis_grad(j, quad.point(q)));
//             pu_grads = pu.node_matrix * pu_grads;
//             real_pu_grads = pu_grads * fv_data.inverse_jacobians[q];
            
            
//             pu_grads.print(cout,"pu_grads");
//             real_pu_grads.print(cout,"real_pu_grads");
            
            //fill regular shape functions
            for (unsigned int k=0; k<dim+1; k++)
            {
//                 unit_grad.zeros();
//                 for (unsigned int l=0; l<dim+1; l++)
//                     unit_grad += fe->basis_grad_vector(l, quad.point(q)) * fe->get_node_matrix()(k,l);
                
//                 unit_grad.print(cout,"unit grad");
                
//                 divs[k] = arma::trace(unit_grad) / fv_data.determinants[q];
//                 DBGCOUT(<< "div=" << fe->basis_div(k, quad.point(q)) << "  |J|=" << fv_data.determinants[q] << "\n");
                divs[k] = fe->basis_div(k, quad.point(q)) / fv_data.determinants[q];
            }
            
            
            //fill enriched shape functions
            unsigned int j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
//                 DBGMSG("interpolant w[%d] q[%d]\n",w,q);
                //compute interpolant
                arma::vec::fixed<spacedim> interpolant; interpolant.zeros();
                double interpolant_div = 0;
                for (unsigned int k=0; k < n_regular_dofs_; k++) {
                    interpolant += fv_data.shape_vectors[q][k] * enr_dof_val[w][k];
                    interpolant_div += divs[k] * enr_dof_val[w][k];
                }
                
//                 quad.real_point(q).print(cout);
//                 DBGMSG("enriched shape value w[%d] q[%d]\n",w,q);
                for (unsigned int k=0; k < pu.n_dofs(); k++)
                {
//                     divs[j] = arma::dot(real_pu_grads.row(k),(enr[w]->vector(fv_data.points[q]) -  interpolant))
//                             + pu_values[k] * (0 - interpolant_div);
//                     divs[j] = arma::dot(real_pu_grads.row(k),enr[w]->vector(fv_data.points[q]));
                    divs[j] = 0;
                    j++;
                }
            }
            
            fv_data.shape_divergence[q] = divs;
        }
    }
}

#endif // FE_RT0_XFEM_HH_