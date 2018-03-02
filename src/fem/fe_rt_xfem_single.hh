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

#ifndef FE_RT0_XFEM_S_HH_
#define FE_RT0_XFEM_S_HH_

#include "fem/finite_element_enriched.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_p.hh"

#include "fem/global_enrichment_func.hh"
#include "fem/fe_values.hh"
#include "quadrature/qxfem.hh"
#include <quadrature/qxfem_factory.hh>
#include "mesh/ref_element.hh"

#include "system/logger.hh"

template <unsigned int dim, unsigned int spacedim> class FE_RT0_XFEM_S;


/**
 * @brief Raviart-Thomas element of order 0.
 *
 * The lowest order Raviart-Thomas finite element with linear basis functions
 * and continuous normal components across element sides.
 */
template <unsigned int dim, unsigned int spacedim>
class FE_RT0_XFEM_S : public FiniteElementEnriched<dim,spacedim>
{
    using FiniteElementEnriched<dim,spacedim>::fe;
    using FiniteElementEnriched<dim,spacedim>::pu;
    using FiniteElementEnriched<dim,spacedim>::enr;
    using FiniteElementEnriched<dim,spacedim>::n_regular_dofs_;
    using FiniteElementEnriched<dim,spacedim>::number_of_dofs;
    
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
public:

    /**
     * @brief Constructor.
     */
    FE_RT0_XFEM_S(FE_RT0<dim,spacedim>* fe,std::vector<EnrichmentPtr> enr);


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
    ~FE_RT0_XFEM_S() {};

private:
    
    void setup_dofs();
    
    /// Awful HACK function for getting quadrature based on Singularity0D or Singularity1D
    std::shared_ptr<QXFEM<dim,spacedim>> qxfem_side(ElementFullIter ele, unsigned int sid);
    
    /// Element data cache for SGFEM interpolation.
    /** Keeps integrals of enrichment functions over sides.
     * This optimization is necessary when outputing on refined output meshes.
     */
    struct EleCache{
        typedef std::vector<std::vector<double>> EnrDofValues;
        std::map<unsigned int, EnrDofValues> enr_dof_values;
    };
    static EleCache ele_cache_;
    
};


template <unsigned int dim, unsigned int spacedim>
typename FE_RT0_XFEM_S<dim,spacedim>::EleCache FE_RT0_XFEM_S<dim,spacedim>::ele_cache_;



template <unsigned int dim, unsigned int spacedim>
FE_RT0_XFEM_S<dim,spacedim>::FE_RT0_XFEM_S(FE_RT0<dim,spacedim>* fe,
                                       std::vector<EnrichmentPtr> enr)
: FiniteElementEnriched<dim,spacedim>(fe,enr)
{
    this->init(false, FEType::FEVectorPiola);
    this->function_space_ = new RT0_space(dim);
    
    setup_dofs();
    // regular + enriched from every singularity
    number_of_dofs = n_regular_dofs_ + enr.size();
    
    this->component_indices_.clear();
    this->nonzero_components_.resize(number_of_dofs, std::vector<bool>(spacedim, true));
}

template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM_S<dim,spacedim>::setup_dofs()
{
    this->dofs_ = fe->dofs_;
    // now add enriched dofs:
    for (unsigned int w=0; w<enr.size(); w++){
        this->dofs_.push_back(Dof(dim, 0,
                                  enr[w]->geometry().dist_vector(arma::vec({0,0,0})),
                                  arma::vec({enr[w]->geometry().effective_surface()}),
                                  DofType::Singular));
    }
}

template <unsigned int dim, unsigned int spacedim>
inline UpdateFlags FE_RT0_XFEM_S<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;

    if (flags & update_values)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    if (flags & update_gradients)
        f |= update_jacobians | update_inverse_jacobians | update_volume_elements | update_quadrature_points;

    return f;
}


template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM_S<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEInternalData &data,
        FEValuesData<dim,spacedim> &fv_data)
{
    fill_fe_values(quad,fv_data);
}

template <unsigned int dim, unsigned int spacedim>
inline void FE_RT0_XFEM_S<dim,spacedim>::fill_fe_values(
        const Quadrature<dim> &quad,
        FEValuesData<dim,spacedim> &fv_data)
{
    ElementFullIter ele = *fv_data.present_cell;
    typedef typename Space<spacedim>::Point Point;
    unsigned int j;

    vector<arma::vec::fixed<spacedim>> normals;
    vector<vector<double>> enr_dof_val;
    
    if (fv_data.update_flags & (update_values | update_divergence | update_gradients))
    {
        normals.resize(RefElement<dim>::n_sides);
//         DBGCOUT("normals\n");
        // compute normals - TODO: this should do mapping, but it does it for fe_side_values..
        for (unsigned int j = 0; j < RefElement<dim>::n_sides; j++){
            normals[j] = arma::trans(fv_data.inverse_jacobians[0])*RefElement<dim>::normal_vector(j);
            normals[j] = normals[j]/norm(normals[j],2);
//             normals[j].print(cout, "internal normal");
        }
        
        auto search = ele_cache_.enr_dof_values.find(ele->index());
        if(search != ele_cache_.enr_dof_values.end()){ // cached
            enr_dof_val = search->second;
//             DBGCOUT(<<"use cached enr dofs on ele " << search->first << "\n");
        }
        else{
//             DBGCOUT(<<"create enr dofs on ele " << ele->index() << "\n");
            enr_dof_val.resize(enr.size());
            for (unsigned int w=0; w<enr.size(); w++)
                enr_dof_val[w].resize(RefElement<dim>::n_sides);
            
            // for SGFEM we need to compute fluxes of glob. enr. function over faces
            for (unsigned int j=0; j<RefElement<dim>::n_sides; j++){
//                 DBGCOUT(<<"qxfem_side on ele " << ele->index() << " s" << j << "\n");
                std::shared_ptr<QXFEM<dim,3>> q_side = qxfem_side(ele, j);
                
                for (unsigned int w=0; w<enr.size(); w++){
                    double val = 0;
//                     DBGCOUT("edge quad [" << j << "]\n");
                    for(unsigned int q=0; q < q_side->size(); q++){
                        val += arma::dot(enr[w]->vector(q_side->real_point(q)),normals[j])
                            // this makes JxW on the triangle side:
//                             * q_side->weight(q)
//                             * side_measure * size_of_ref_element;
                               * q_side->JxW(q);
                    }
                    enr_dof_val[w][j] = val;
//                     DBGVAR(val);
//                     cout << setprecision(16) << "val  " << val << endl;
                }
            }
            //copy to cache
            ele_cache_.enr_dof_values[ele->index()] = enr_dof_val;
        }
//         // for SGFEM we need to compute fluxes of glob. enr. function over faces
//         ASSERT(dim == 2);
//         const uint qsize=100;
//         QMidpoint qside(qsize);
//         MappingP1<2,3> map;
//         arma::mat ele_mat = map.element_map(*ele);
//         
//         for (unsigned int w=0; w<enr.size(); w++){
//             enr_dof_val[w].resize(RefElement<dim>::n_sides);
//             for (unsigned int j=0; j<RefElement<dim>::n_sides; j++){
//                 
//                 Quadrature<2> quad_side(qside, j, *ele->permutation_idx_);
//                     
//                 double side_measure = ele->side(j)->measure();
//                 double val = 0;
// //                 DBGCOUT("edge quad [" << j << "]\n");
//                 for(unsigned int q=0; q < quad_side.size(); q++){
//                     
//                     arma::vec3 real_point = map.project_unit_to_real(RefElement<2>::local_to_bary(quad_side.point(q)), ele_mat);
// //                     arma::vec3 real_point = ele_mat.cols(1,dim) * quad_side.point(q) + ele_mat.col(0);
// //                     cout << real_point(0) << " " << real_point(1) << " " << real_point(2) << "\n";
//                     val += arma::dot(enr[w]->vector(real_point),normals[j])
//                         // this makes JxW on the triangle side:
//                         * quad_side.weight(q)
//                         * side_measure;
//                 }
//                 enr_dof_val[w][j] = val;
// //                 DBGVAR(val);
// //                 cout << setprecision(16) << "val  " << val << endl;
//             }
//         }
      
    }
    
    
    // shape values
    if (fv_data.update_flags & update_values)
    {
        vector<arma::vec::fixed<spacedim> > vectors(number_of_dofs);
        
        for (unsigned int q = 0; q < quad.size(); q++)
        {
            //fill regular shape functions
            arma::mat::fixed<dim+1,dim> raw_values;
            arma::mat::fixed<dim+1,dim> shape_values;
            
//             DBGCOUT("regular shape vals\n");
            for (j=0; j<n_regular_dofs_; j++)
                for (unsigned int c=0; c<dim; c++)
                    raw_values(j,c) = fe->basis_value(j, quad.point(q), c);
            
//             DBGCOUT("times node matrix\n");
            shape_values = fe->get_node_matrix() * raw_values;
            
//             DBGCOUT("times jacobians\n");
            for (j=0; j<n_regular_dofs_; j++)
                vectors[j] = fv_data.jacobians[q] * arma::trans(shape_values.row(j)) / fv_data.determinants[q];
            
            //fill enriched shape functions
//             DBGCOUT("enriched shape vals\n");
            j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
//                 DBGMSG("interpolant w[%d] q[%d]\n",w,q);
                //compute interpolant
                arma::vec::fixed<spacedim> interpolant;
                interpolant.zeros();
                for (unsigned int k=0; k < n_regular_dofs_; k++)
                    interpolant += vectors[k] * enr_dof_val[w][k];
                
                vectors[j] =  enr[w]->vector(fv_data.points[q]) - interpolant;
                j++;
            }   
            
//             DBGCOUT("enriched shape vals - decomp on components\n");
            for (unsigned int k=0; k<number_of_dofs; k++)
                for (unsigned int c=0; c<dim; c++)
                    fv_data.shape_values[q][k*spacedim+c] = vectors[k][c];
        }
    }

    // divergence
//     DBGCOUT("divergence\n");
    if (fv_data.update_flags & update_divergence)
    {
        vector<double> divs(number_of_dofs);
        
        arma::mat::fixed<dim,dim> unit_grad;
        arma::mat::fixed<spacedim,spacedim> real_grad;
        for (unsigned int q = 0; q < quad.size(); q++)
        {   
            //fill regular shape functions
//             DBGCOUT("regular shape divs\n");
            for (unsigned int k=0; k<dim+1; k++)
            {
                // grads on ref element
                unit_grad.zeros();
                  
                for (unsigned int l=0; l<dim+1; l++)
                    for (unsigned int c=0; c<dim; c++)
                        unit_grad.col(c) += fe->basis_grad(l, quad.point(q),c) * fe->get_node_matrix()(k,l);
                
                // map grads on real element
                real_grad = fv_data.jacobians[q] * unit_grad * fv_data.inverse_jacobians[q]/fv_data.determinants[q];

                // compute div as a trace
                divs[k] = arma::trace(real_grad);
//                 DBGCOUT(<< "div=" << divs[k] << "\n");
            }
            
            
            //fill enriched shape functions
//             DBGCOUT("enriched shape divs\n");
            j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
                //compute interpolant
                double interpolant_div = 0;
                for (unsigned int k=0; k < n_regular_dofs_; k++) {
                    interpolant_div += divs[k] * enr_dof_val[w][k];
                }
                
                divs[j] = - interpolant_div;
                j++;
            }
            
            fv_data.shape_divergence[q] = divs;
        }
    }
}

template<>
inline std::shared_ptr<QXFEM<1,3>> FE_RT0_XFEM_S<1,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    return std::make_shared<QXFEM<1,3>>();
}

template<>
inline std::shared_ptr<QXFEM<2,3>> FE_RT0_XFEM_S<2,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    std::vector<std::shared_ptr<Singularity<0>>> vec(enr.size());
    for(unsigned int w=0; w<enr.size(); w++){
        vec[w] = static_pointer_cast<Singularity<0>>(enr[w]);
    }
    QXFEMFactory qfact(12);
    return qfact.create_side_singular(vec, ele, sid);
}

template<>
inline std::shared_ptr<QXFEM<3,3>> FE_RT0_XFEM_S<3,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    std::vector<std::shared_ptr<Singularity<1>>> vec(enr.size());
    for(unsigned int w=0; w<enr.size(); w++){
        vec[w] = static_pointer_cast<Singularity<1>>(enr[w]);
    }
    QXFEMFactory qfact(12);
    return qfact.create_side_singular(vec, ele, sid);
}

#endif // FE_RT0_XFEM_S_HH_
