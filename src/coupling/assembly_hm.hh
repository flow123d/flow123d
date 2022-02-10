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
 * @file    assembly_elasticity.hh
 * @brief
 */

#ifndef ASSEMBLY_HM_HH_
#define ASSEMBLY_HM_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "coupling/hm_iterative.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class FlowPotentialAssemblyHM : public AssemblyBase<dim>
{
public:
    typedef typename HM_Iterative::EqFields EqFields;
    typedef void EqData;

    static constexpr const char * name() { return "FlowPotentialAssemblyHM"; }

    /// Constructor.
    FlowPotentialAssemblyHM(EqFields *eq_fields, EqData *)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields) {
        this->active_integrals_ = (ActiveIntegrals::bulk);
        this->used_fields_ += eq_fields_->alpha;
        this->used_fields_ += eq_fields_->density;
        this->used_fields_ += eq_fields_->gravity;
    }

    /// Destructor.
    ~FlowPotentialAssemblyHM() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_low_ = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        fe_values_.initialize(*this->quad_, *fe_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);
        fe_values_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
        fe_values_sub_.initialize(*this->quad_low_, *fe_low_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);

        n_dofs_ = fe_->n_dofs();
        n_dofs_sub_ = fe_low_->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2);
        side_dof_indices_[0].resize(n_dofs_sub_);  // index 0 = element with lower dimension,
        side_dof_indices_[1].resize(n_dofs_);      // index 1 = side of element with higher dimension
        local_matrix_.resize(n_dofs_*n_dofs_);
        local_matrix_ngh_.resize(2);
        for (uint m=0; m<2; ++m) {
            local_matrix_ngh_[m].resize(2);
            for (uint n=0; n<2; ++n)
                local_matrix_ngh_[m][n].resize(n_dofs_*n_dofs_);
        }
        vec_view_ = &fe_values_.vector_view(0);
        vec_view_side_ = &fe_values_side_.vector_view(0);
        if (dim>1) vec_view_sub_ = &fe_values_sub_.vector_view(0);
    }


    /// Assemble integral over element
    // inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    // {
    //     if (cell.dim() != dim) return;

    //     ElementAccessor<3> elm_acc = cell.elm();

    //     fe_values_.reinit(elm_acc);
    //     cell.get_dof_indices(dof_indices_);

    //     // assemble the local stiffness matrix
    //     for (unsigned int i=0; i<n_dofs_; i++)
    //         for (unsigned int j=0; j<n_dofs_; j++)
    //             local_matrix_[i*n_dofs_+j] = 0;

    //     unsigned int k=0;
    //     for (auto p : this->bulk_points(element_patch_idx) )
    //     {
        //     for (unsigned int i=0; i<n_dofs_; i++)
        //     {
        //         for (unsigned int j=0; j<n_dofs_; j++)
        //             local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
        //                                         2*eq_fields_->lame_mu(p)*arma::dot(vec_view_->sym_grad(j,k), vec_view_->sym_grad(i,k))
        //                                         + eq_fields_->lame_lambda(p)*vec_view_->divergence(j,k)*vec_view_->divergence(i,k)
        //                                        )*fe_values_.JxW(k);
        //     }
        //     k++;
        // }
        // eq_data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    // }


private:
    inline arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
    {
      arma::mat33 mt = m - m*arma::kron(n,n.t());
      return mt;
    }



    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)
    FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fe_values_side_;                              ///< FEValues of side object
    FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    vector<vector<vector<PetscScalar>>> local_matrix_ngh_;    ///< Auxiliary vectors for assemble ngh integral
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in boundary / neighbour calculation.
    const FEValuesViews::Vector<3> * vec_view_sub_;           ///< Vector view of low dim element in neighbour calculation.

    template <template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};




#endif /* ASSEMBLY_HM_HH_ */

