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
 * @file    fe_values.cc
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel
 */

#include "fem/mapping.hh"
#include "quadrature/quadrature.hh"
#include <quadrature/qxfem.hh>
#include "fem/finite_element.hh"
#include "fem/xfe_values.hh"

#include "fem/singularity.hh"
#include "fem/xfem_element_data.hh"
#include "quadrature/qxfem_factory.hh"


using namespace arma;
using namespace std;

// template <unsigned int dim, unsigned int spacedim>
// void XFEValues<dim,spacedim>::setup_dofs()
// {
//     this->dofs_ = fe->dofs_;
//     // now add enriched dofs:
//     for (unsigned int w=0; w<enr.size(); w++){
//         this->dofs_.push_back(Dof(dim, 0,
//                                   enr[w]->geometry().dist_vector(arma::vec({0,0,0})),
//                                   arma::vec({enr[w]->geometry().effective_surface()}),
//                                   DofType::Singular));
//     }
// }

// template <unsigned int dim, unsigned int spacedim>
// inline void FE_RT0_XFEM_S<dim,spacedim>::fill_fe_values(
//         const Quadrature<dim> &quad,
//         FEInternalData &data,
//         FEValuesData<dim,spacedim> &fv_data)
// {
//     fill_fe_values(quad,fv_data);
// }

template <unsigned int dim, unsigned int spacedim>
typename XFEValues<dim,spacedim>::EleCache XFEValues<dim,spacedim>::ele_cache;

template<unsigned int dim,unsigned int spacedim>
XFEValues<dim,spacedim>::XFEValues(
            Mapping<dim,spacedim> &mapping,
            FiniteElement<dim> &fe,
            FiniteElement<dim> &pu,
            UpdateFlags flags)
: FEValuesBase<dim, spacedim>(),
pu(&pu),
n_regular_dofs_ (fe.n_dofs()),
n_enriched_dofs_(0)
{
    this->mapping = &mapping;
    this->fe = &fe;
    this->data.update_flags = flags;
}

template<unsigned int dim,unsigned int spacedim>
void XFEValues<dim,spacedim>::reinit(ElementFullIter & ele,
                                     XFEMElementData<dim,spacedim> & xdata,
                                     Quadrature<dim> &_quadrature)
{
    ASSERT_EQ_DBG( dim, ele->dim() );
    this->data.present_cell = &ele;
    this->quadrature = &_quadrature;
    
    enr = xdata.enrichment_func_vec();
    n_enriched_dofs_ = enr.size()*pu->n_dofs();
    
    this->allocate(*this->mapping, *this->quadrature, *this->fe, n_regular_dofs_ + n_enriched_dofs_, this->data.update_flags);
    
    if(typeid(_quadrature) == typeid(QXFEM<dim,spacedim>)){
        QXFEM<dim,spacedim>* q = static_cast<QXFEM<dim,spacedim>*>(&_quadrature);
        
        this->data.update_flags = this->data.update_flags & (~update_quadrature_points);
        //TODO: think of way to avoid copying the whole vector
        this->data.points = q->get_real_points();
    }
    
    // precompute the maping data and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_values(ele,
                            *this->quadrature,
                            *this->mapping_data,
                            this->data);

    this->fill_data(*this->fe_data);
}

template<unsigned int dim, unsigned int spacedim>
void XFEValues<dim,spacedim>::fill_data(const FEInternalData &fe_data)
{
//     DebugOut() << "XFEValues::fill_data\n";
    switch (fe->type_) {
        case FEScalar:
            fill_scalar_xfem_single();
            break;
//         case FEVectorContravariant:
//             fill_vec_contravariant_data(fe_data);
//             break;
        case FEVectorPiola:
            fill_vec_piola_xfem_single();
            break;
//         case FEMixedSystem:
//             fill_system_data(fe_data);
//             break;
        default:
            ASSERT(false).error("Not implemented.");
    }
}

template <unsigned int dim, unsigned int spacedim>
void XFEValues<dim,spacedim>::fill_scalar_xfem_single()
{
//     DebugOut() << "XFEValues::fill_scalar_xfem_single\n";
    ASSERT_DBG(fe->type_ == FEScalar);
    
    // shape values
    if (this->data.update_flags & update_values)
        for (unsigned int i = 0; i < this->quadrature->size(); i++)
            for (unsigned int j = 0; j < fe->n_dofs(); j++)
                this->data.shape_values[i][j] = fe->shape_value(j, this->quadrature->point(i), 0);
        
    // shape gradients
    arma::mat grad(dim, 0);
    if (this->data.update_flags & update_gradients)
        for (unsigned int i = 0; i < this->quadrature->size(); i++)
            for (unsigned int j = 0; j < fe->n_dofs(); j++)
            {
                grad.zeros();
                for (unsigned int c=0; c<fe->n_components(); c++)
                    grad.col(c) += fe->shape_grad(j, this->quadrature->point(i), c);
        
                this->data.shape_gradients[i][j] = trans(this->data.inverse_jacobians[i]) * grad;
            }
}

template <unsigned int dim, unsigned int spacedim>
void XFEValues<dim,spacedim>::fill_vec_piola_xfem_single()
{
//     DebugOut() << "XFEValues::fill_vec_piola_xfem_single\n";
    ASSERT_DBG(fe->type_ == FEVectorPiola);
    
    ElementFullIter ele = *this->data.present_cell;
    typedef typename Space<spacedim>::Point Point;
    unsigned int j;

    vector<arma::vec::fixed<spacedim>> normals;
    vector<vector<double>> enr_dof_val;
    
    if (this->data.update_flags & (update_values | update_divergence | update_gradients))
    {
        ASSERT_DBG(this->data.inverse_jacobians.size() > 0);
        
        normals.resize(RefElement<dim>::n_sides);
//         DBGCOUT("normals\n");
        // compute normals - TODO: this should do mapping, but it does it for fe_side_values..
        for (unsigned int j = 0; j < RefElement<dim>::n_sides; j++){
            normals[j] = arma::trans(this->data.inverse_jacobians[0])*RefElement<dim>::normal_vector(j);
            normals[j] = normals[j]/norm(normals[j],2);
//             normals[j].print(cout, "internal normal");
        }
        
        auto search = ele_cache.enr_dof_values.find(ele->index());
        if(search != ele_cache.enr_dof_values.end()){ // cached
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
            ele_cache.enr_dof_values[ele->index()] = enr_dof_val;
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
    DBGCOUT("shape values\n");
    if (this->data.update_flags & update_values)
    {
        ASSERT_DBG(this->data.jacobians.size() == this->quadrature->size());
        ASSERT_DBG(this->data.determinants.size() == this->quadrature->size());
        ASSERT_DBG(this->data.points.size() == this->quadrature->size());
        
        vector<arma::vec::fixed<spacedim> > vectors(this->n_dofs_);
        
        for (unsigned int q = 0; q < this->n_points(); q++)
        {
            //fill regular shape functions
//             DBGCOUT("regular shape vals\n");
            arma::mat shape_values(this->fe->n_dofs(), this->fe->n_components());
            for (unsigned int j=0; j<this->fe->n_dofs(); j++)
                for (unsigned int c=0; c<this->fe->n_components(); c++)
                    shape_values(j,c) = this->fe->shape_value(j, this->quadrature->point(q), c);
            
//             DBGCOUT("times jacobians\n");
            for (j=0; j<n_regular_dofs_; j++)
                vectors[j] = this->data.jacobians[q] * arma::trans(shape_values.row(j)) / this->data.determinants[q];
            
            //fill enriched shape functions
//             DBGCOUT("enriched shape vals\n");
            j = n_regular_dofs_;
            for (unsigned int w=0; w<enr.size(); w++)
            {
                //compute interpolant
                arma::vec::fixed<spacedim> interpolant;
                interpolant.zeros();
                for (unsigned int k=0; k < n_regular_dofs_; k++)
                    interpolant += vectors[k] * enr_dof_val[w][k];
                
                vectors[j] =  enr[w]->vector(this->data.points[q]) - interpolant;
                j++;
            }   
            
//             DBGCOUT("enriched shape vals - decomp on components\n");
            for (unsigned int k=0; k<this->n_dofs_; k++)
                for (unsigned int c=0; c<spacedim; c++)
                    this->data.shape_values[q][k*spacedim+c] = vectors[k][c];
        }
    }

    // divergence
    DBGCOUT("divergence\n");
    if (this->data.update_flags & update_divergence)
    {
        ASSERT_DBG(this->data.jacobians.size() == this->quadrature->size());
        ASSERT_DBG(this->data.inverse_jacobians.size() == this->quadrature->size());
        ASSERT_DBG(this->data.determinants.size() == this->quadrature->size());
        ASSERT_DBG(this->data.points.size() == this->quadrature->size());
        
        vector<double> divs(this->n_dofs_);
        
        arma::mat::fixed<dim,dim> unit_grad;
        arma::mat::fixed<spacedim,spacedim> real_grad;
        for (unsigned int q = 0; q < this->n_points(); q++)
        {   
            //fill regular shape functions
//             DBGCOUT("regular shape divs\n");
            for (unsigned int k=0; k<dim+1; k++)
            {
//                 // grads on ref element
                unit_grad.zeros();
                for (unsigned int c=0; c<fe->n_components(); c++)
                    unit_grad.col(c) += this->fe->shape_grad(k, this->quadrature->point(q),c);
                
                // map grads on real element
                real_grad = trans(this->data.inverse_jacobians[q]) * unit_grad * trans(this->data.jacobians[q])/this->data.determinants[q];

//                 for (unsigned int c=0; c<spacedim; c++)
//                     data.shape_gradients[i][j*spacedim+c] = grads.col(c);
                
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
            
            this->data.shape_divergence[q] = divs;
        }
    }
}

template<>
inline std::shared_ptr<QXFEM<1,3>> XFEValues<1,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    return std::make_shared<QXFEM<1,3>>();
}

template<>
inline std::shared_ptr<QXFEM<2,3>> XFEValues<2,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    std::vector<std::shared_ptr<Singularity<0>>> vec(enr.size());
    for(unsigned int w=0; w<enr.size(); w++){
        vec[w] = static_pointer_cast<Singularity<0>>(enr[w]);
    }
    QXFEMFactory qfact(12);
    return qfact.create_side_singular(vec, ele, sid);
}

template<>
inline std::shared_ptr<QXFEM<3,3>> XFEValues<3,3>::qxfem_side(ElementFullIter ele, unsigned int sid)
{
    std::vector<std::shared_ptr<Singularity<1>>> vec(enr.size());
    for(unsigned int w=0; w<enr.size(); w++){
        vec[w] = static_pointer_cast<Singularity<1>>(enr[w]);
    }
    QXFEMFactory qfact(12);
    return qfact.create_side_singular(vec, ele, sid);
}



template class XFEValues<1,3>;
template class XFEValues<2,3>;
template class XFEValues<3,3>;

// template class FESideValues<1,3>;
// template class FESideValues<2,3>;
// template class FESideValues<3,3>;

