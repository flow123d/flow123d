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

#ifndef ASSEMBLY_ELASTICITY_HH_
#define ASSEMBLY_ELASTICITY_HH_

#include <algorithm>
#include <cmath>
#include <map>
#include <queue>

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "mechanics/elasticity.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"
#include "la/linsys_PERMON.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class StiffnessAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssemblyElasticity"; }

    /// Constructor.
    StiffnessAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->lame_mu;
        this->used_fields_ += eq_fields_->lame_lambda;
        this->used_fields_ += eq_fields_->dirichlet_penalty;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->fracture_sigma;
    }

    /// Destructor.
    ~StiffnessAssemblyElasticity() {}

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
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        dof_indices_ = cell.get_loc_dof_indices();

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        unsigned int k=0;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
                                                2*eq_fields_->lame_mu(p)*arma::dot(vec_view_->sym_grad(j,k), vec_view_->sym_grad(i,k))
                                                + eq_fields_->lame_lambda(p)*vec_view_->divergence(j,k)*vec_view_->divergence(i,k)
                                               )*fe_values_.JxW(k);
            }
            k++;
        }
        eq_data_->ls->mat_set_values_local(n_dofs_, dof_indices_.memptr(), n_dofs_, dof_indices_.memptr(), &(local_matrix_[0]));
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (eq_data_->dirichlet_by_eq) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dof_indices_ = dh_cell.get_loc_dof_indices();
        fe_values_side_.reinit(side);

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqFields::bc_type_displacement)
        {
            unsigned int k=0;
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k),vec_view_side_->value(j,k)) * fe_values_side_.JxW(k);
                k++;
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            unsigned int k=0;
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k), fe_values_side_.normal_vector(k)) *
                                arma::dot(vec_view_side_->value(j,k), fe_values_side_.normal_vector(k)) * fe_values_side_.JxW(k);
                k++;
            }
        }

        eq_data_->ls->mat_set_values_local(n_dofs_, dof_indices_.memptr(), n_dofs_, dof_indices_.memptr(), &(local_matrix_[0]));
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        if (!cell_lower_dim.is_own()) return;

		side_dof_indices_[0] = cell_lower_dim.get_loc_dof_indices();
		fe_values_sub_.reinit(cell_lower_dim.elm());
        side_dof_indices_[1] = cell_higher_dim.get_loc_dof_indices();
		fe_values_side_.reinit(neighb_side.side());

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int m=0; m<2; ++m)
                for (unsigned int i=0; i<n_dofs_*n_dofs_; i++)
                    local_matrix_ngh_[n][m][i] = 0;

        // set transmission conditions
        unsigned int k=0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = fe_values_side_.normal_vector(k);

            for (int n=0; n<2; n++)
            {
                for (unsigned int i=0; i<n_dofs_ngh_[n]; i++)
                {
                    arma::vec3 vi = (n==0) ? arma::zeros(3) : vec_view_side_->value(i,k);
                    arma::vec3 vf = (n==1) ? arma::zeros(3) : vec_view_sub_->value(i,k);
                    arma::mat33 gvft = (n==0) ? mat_t(vec_view_sub_->grad(i,k),nv) : arma::zeros(3,3);
                    arma::mat33 gvi = (n==0) ? arma::zeros(3,3) : vec_view_side_->grad(i,k);
                    if (eq_data_->fix_nullspace) {
                        vi = vi - eq_fields_->cross_section(p_low)/2 * (gvi*nv);
                    }
                    double divvft = (n==0) ? arma::trace(gvft) : 0;

                    for (int m=0; m<2; m++)
                    {
                        for (unsigned int j=0; j<n_dofs_ngh_[m]; j++) {
                            arma::vec3 ui = (m==0) ? arma::zeros(3) : vec_view_side_->value(j,k);
                            arma::vec3 uf = (m==1) ? arma::zeros(3) : vec_view_sub_->value(j,k);
                            arma::mat33 guft = (m==0) ? mat_t(vec_view_sub_->grad(j,k),nv) : arma::zeros(3,3);
                            arma::mat33 gui = (m==0) ? arma::zeros(3,3) : vec_view_side_->grad(j,k);
                            double divuft = (m==0) ? arma::trace(guft) : 0;
                            if (eq_data_->fix_nullspace) {
                                ui = ui - eq_fields_->cross_section(p_low)/2 * (gui*nv);
                            }

                            local_matrix_ngh_[n][m][i*n_dofs_ngh_[m] + j] +=
                                    eq_fields_->fracture_sigma(p_low)*(
                                        2/eq_fields_->cross_section(p_low)*(
                                            eq_fields_->lame_mu(p_low)*arma::dot(vf-vi,uf-ui)
                                          +(eq_fields_->lame_mu(p_low)+eq_fields_->lame_lambda(p_low))*arma::dot(uf-ui,nv)*arma::dot(vf-vi,nv)
                                        )
                                        + eq_fields_->lame_mu(p_low)*( arma::dot(vf-vi,guft.t()*nv) + arma::dot(uf-ui,gvft.t()*nv) )
                                        + eq_fields_->lame_lambda(p_low)*( divuft*arma::dot(vf-vi,nv) + divvft*arma::dot(uf-ui,nv) )
                                    )*fe_values_sub_.JxW(k) * 2/cell_lower_dim.elm()->n_neighs_vb();
                            
                            // artificial tangential stiffness on adjacent dofs
                            if (eq_data_->fix_nullspace) {
                                local_matrix_ngh_[n][m][i*n_dofs_ngh_[m] + j] +=
                                    eq_fields_->fracture_sigma(p_low)*(
                                        eq_data_->fix_nullspace ? ( eq_fields_->cross_section(p_low)*(
                                                2*eq_fields_->lame_mu(p_low)*arma::dot((gui+gui.t())/2, (gvi+gvi.t())/2)
                                                + eq_fields_->lame_lambda(p_low)*arma::trace(gvi)*arma::trace(gui)
                                               )
                                            ):0

                                    )*fe_values_sub_.JxW(k) * 2/cell_lower_dim.elm()->n_neighs_vb();
                            }
                        }
                    }

                }
            }
        	k++;
        }

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int m=0; m<2; ++m)
                eq_data_->ls->mat_set_values_local(n_dofs_ngh_[n], side_dof_indices_[n].memptr(), n_dofs_ngh_[m], side_dof_indices_[m].memptr(), &(local_matrix_ngh_[n][m][0]));
    }



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
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)
    FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fe_values_side_;                              ///< FEValues of side object
    FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

    LocDofVec dof_indices_;                             ///< Vector of global DOF indices
    vector<LocDofVec > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    vector<vector<vector<PetscScalar>>> local_matrix_ngh_;    ///< Auxiliary vectors for assemble ngh integral
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in boundary / neighbour calculation.
    const FEValuesViews::Vector<3> * vec_view_sub_;           ///< Vector view of low dim element in neighbour calculation.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class RhsAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "RhsAssemblyElasticity"; }

    /// Constructor.
    RhsAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->load;
        this->used_fields_ += eq_fields_->potential_load;
        this->used_fields_ += eq_fields_->ref_potential_load;
        this->used_fields_ += eq_fields_->fracture_sigma;
        this->used_fields_ += eq_fields_->dirichlet_penalty;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_displacement;
        this->used_fields_ += eq_fields_->bc_traction;
        this->used_fields_ += eq_fields_->bc_stress;
        this->used_fields_ += eq_fields_->initial_stress;
    }

    /// Destructor.
    ~RhsAssemblyElasticity() {}

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
        fe_values_bdr_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
        fe_values_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_normal_vectors);
        fe_values_sub_.initialize(*this->quad_low_, *fe_low_,
                update_values | update_JxW_values | update_quadrature_points);
        n_dofs_ = fe_->n_dofs();
        n_dofs_sub_ = fe_low_->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2);
        side_dof_indices_[0].resize(n_dofs_sub_);  // index 0 = element with lower dimension,
        side_dof_indices_[1].resize(n_dofs_);      // index 1 = side of element with higher dimension
        local_rhs_.resize(n_dofs_);
        local_rhs_ngh_.resize(2);
        for (uint n=0; n<2; ++n) local_rhs_ngh_[n].resize(n_dofs_);
        vec_view_ = &fe_values_.vector_view(0);
        vec_view_bdr_ = &fe_values_bdr_side_.vector_view(0);
        vec_view_side_ = &fe_values_side_.vector_view(0);
        if (dim>1) vec_view_sub_ = &fe_values_sub_.vector_view(0);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        //local_source_balance_vector.assign(n_dofs_, 0);
        //local_source_balance_rhs.assign(n_dofs_, 0);

        // compute sources
        unsigned int k=0;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += (
                                 arma::dot(eq_fields_->load(p), vec_view_->value(i,k))
                                 -eq_fields_->potential_load(p)*vec_view_->divergence(i,k)
                                 -arma::dot(eq_fields_->initial_stress(p), vec_view_->grad(i,k))
                                )*eq_fields_->cross_section(p)*fe_values_.JxW(k);
            ++k;
        }
        eq_data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));

//         for (unsigned int i=0; i<n_dofs_; i++)
//         {
//             for (unsigned int k=0; k<qsize_; k++) // point range
//                 local_source_balance_vector[i] -= 0;//sources_sigma[k]*fe_values_[vec_view_].value(i,k)*fe_values_.JxW(k);
//
//             local_source_balance_rhs[i] += local_rhs_[i];
//         }
//         balance_->add_source_matrix_values(subst_idx, elm_acc.region().bulk_idx(), dof_indices_, local_source_balance_vector);
//         balance_->add_source_vec_values(subst_idx, elm_acc.region().bulk_idx(), dof_indices_, local_source_balance_rhs);
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);
        fe_values_bdr_side_.reinit(side);

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        // local_flux_balance_vector.assign(n_dofs_, 0);
        // local_flux_balance_rhs = 0;

        unsigned int k = 0;

        // addtion from initial stress
        for (auto p : this->boundary_points(cell_side) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += eq_fields_->cross_section(p) *
                        arma::dot(( eq_fields_->initial_stress(p) * fe_values_bdr_side_.normal_vector(k)),
                                    vec_view_bdr_->value(i,k)) *
                        fe_values_bdr_side_.JxW(k);
            ++k;
        }

        k = 0;
        if (bc_type == EqFields::bc_type_displacement && !eq_data_->dirichlet_by_eq)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(eq_fields_->bc_displacement(p_bdr), vec_view_bdr_->value(i,k)) *
					        fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal && !eq_data_->dirichlet_by_eq)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                            arma::dot(eq_fields_->bc_displacement(p_bdr), fe_values_bdr_side_.normal_vector(k)) *
                            arma::dot(vec_view_bdr_->value(i,k), fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_traction)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(vec_view_bdr_->value(i,k), eq_fields_->bc_traction(p_bdr) + eq_fields_->ref_potential_load(p) * fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_stress)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    // stress is multiplied by inward normal to obtain traction
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(vec_view_bdr_->value(i,k), -eq_fields_->bc_stress(p_bdr)*fe_values_bdr_side_.normal_vector(k)
                            + eq_fields_->ref_potential_load(p) * fe_values_bdr_side_.normal_vector(k))
                            * fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        eq_data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));


//             balance_->add_flux_matrix_values(subst_idx, loc_b, side_dof_indices, local_flux_balance_vector);
//             balance_->add_flux_vec_value(subst_idx, loc_b, local_flux_balance_rhs);
		// ++loc_b;
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

		cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
		ElementAccessor<3> cell_sub = cell_lower_dim.elm();
		fe_values_sub_.reinit(cell_sub);

		DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
		cell_higher_dim.get_dof_indices(side_dof_indices_[1]);
		fe_values_side_.reinit(neighb_side.side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_ngh_[n][i] = 0;

        // set transmission conditions
        unsigned int k=0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = fe_values_side_.normal_vector(k);

            for (int n=0; n<2; n++)
            {
                if (!own_element_id[n]) continue;

                for (unsigned int i=0; i<n_dofs_ngh_[n]; i++)
                {
                    arma::vec3 vi = (n==0) ? arma::zeros(3) : vec_view_side_->value(i,k);
                    arma::vec3 vf = (n==1) ? arma::zeros(3) : vec_view_sub_->value(i,k);

                    local_rhs_ngh_[n][i] -= eq_fields_->fracture_sigma(p_low) * eq_fields_->cross_section(p_high) *
                            arma::dot(vf-vi, eq_fields_->potential_load(p_high) * nv) * fe_values_sub_.JxW(k);
                }
            }
            ++k;
        }

        for (unsigned int n=0; n<2; ++n)
            eq_data_->ls->rhs_set_values(n_dofs_ngh_[n], side_dof_indices_[n].data(), &(local_rhs_ngh_[n][0]));
    }



private:
    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)
    FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fe_values_bdr_side_;                          ///< FEValues of side (boundary integral) object
    FEValues<3> fe_values_side_;                              ///< FEValues of side (neighbour integral) object
    FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for assemble methods
    vector<vector<PetscScalar>> local_rhs_ngh_;               ///< Auxiliary vectors for assemble ngh integral
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_bdr_;           ///< Vector view in boundary calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in neighbour calculation.
    const FEValuesViews::Vector<3> * vec_view_sub_;           ///< Vector view of low dim element in neighbour calculation.


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

template <unsigned int dim>
class OutpuFieldsAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "OutpuFieldsAssemblyElasticity"; }

    /// Constructor.
    OutpuFieldsAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->lame_mu;
        this->used_fields_ += eq_fields_->lame_lambda;
        this->used_fields_ += eq_fields_->initial_stress;
    }

    /// Destructor.
    ~OutpuFieldsAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fv_.initialize(*this->quad_, *fe_,
        		update_values | update_gradients | update_quadrature_points);
        fsv_.initialize(*this->quad_low_, *fe_,
        		update_values | update_normal_vectors | update_quadrature_points);
        n_dofs_ = fe_->n_dofs();
        vec_view_ = &fv_.vector_view(0);
        //        if (dim>1) ??
        vec_view_side_ = &fsv_.vector_view(0);

        output_vec_ = eq_fields_->output_field_ptr->vec();
        output_stress_vec_ = eq_fields_->output_stress_ptr->vec();
        output_von_mises_stress_vec_ = eq_fields_->output_von_mises_stress_ptr->vec();
        output_mean_stress_vec_ = eq_fields_->output_mean_stress_ptr->vec();
        output_cross_sec_vec_ = eq_fields_->output_cross_section_ptr->vec();
        output_div_vec_ = eq_fields_->output_div_ptr->vec();
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;
        DHCellAccessor cell_tensor = cell.cell_with_other_dh(eq_data_->dh_tensor_.get());
        DHCellAccessor cell_scalar = cell.cell_with_other_dh(eq_data_->dh_scalar_.get());

        auto elm = cell.elm();

        fv_.reinit(elm);
        dof_indices_        = cell.get_loc_dof_indices();
        dof_indices_scalar_ = cell_scalar.get_loc_dof_indices();
        dof_indices_tensor_ = cell_tensor.get_loc_dof_indices();

        auto p = *( this->bulk_points(element_patch_idx).begin() );

        arma::mat33 stress = eq_fields_->initial_stress(p);
        double div = 0;
        for (unsigned int i=0; i<n_dofs_; i++)
        {
            stress += (2*eq_fields_->lame_mu(p)*vec_view_->sym_grad(i,0) + eq_fields_->lame_lambda(p)*vec_view_->divergence(i,0)*arma::eye(3,3))
                    * output_vec_.get(dof_indices_[i]);
            div += vec_view_->divergence(i,0)*output_vec_.get(dof_indices_[i]);
        }

        arma::mat33 stress_dev = stress - arma::trace(stress)/3*arma::eye(3,3);
        double von_mises_stress = sqrt(1.5*arma::dot(stress_dev, stress_dev));
        double mean_stress = arma::trace(stress) / 3;
        output_div_vec_.add(dof_indices_scalar_[0], div);

        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<3; j++)
                output_stress_vec_.add( dof_indices_tensor_[i*3+j], stress(i,j) );
        output_von_mises_stress_vec_.set( dof_indices_scalar_[0], von_mises_stress );
        output_mean_stress_vec_.set( dof_indices_scalar_[0], mean_stress );

        output_cross_sec_vec_.add( dof_indices_scalar_[0], eq_fields_->cross_section(p) );
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        normal_displacement_ = 0;
        normal_stress_.zeros();

        DHCellAccessor cell_higher_dim = neighb_side.cell();
        DHCellAccessor cell_tensor = cell_lower_dim.cell_with_other_dh(eq_data_->dh_tensor_.get());
        DHCellAccessor cell_scalar = cell_lower_dim.cell_with_other_dh(eq_data_->dh_scalar_.get());
        fsv_.reinit(neighb_side.side());

        dof_indices_ = cell_higher_dim.get_loc_dof_indices();
        auto p_high = *( this->coupling_points(neighb_side).begin() );
        auto p_low = p_high.lower_dim(cell_lower_dim);

        for (unsigned int i=0; i<n_dofs_; i++)
        {
            normal_displacement_ -= arma::dot(vec_view_side_->value(i,0)*output_vec_.get(dof_indices_[i]), fsv_.normal_vector(0));
            arma::mat33 grad = -arma::kron(vec_view_side_->value(i,0)*output_vec_.get(dof_indices_[i]), fsv_.normal_vector(0).t()) / eq_fields_->cross_section(p_low);
            normal_stress_ += eq_fields_->lame_mu(p_low)*(grad+grad.t()) + eq_fields_->lame_lambda(p_low)*arma::trace(grad)*arma::eye(3,3);
        }

        LocDofVec dof_indices_scalar_ = cell_scalar.get_loc_dof_indices();
        LocDofVec dof_indices_tensor_ = cell_tensor.get_loc_dof_indices();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<3; j++)
                output_stress_vec_.add( dof_indices_tensor_[i*3+j], normal_stress_(i,j) );
        output_cross_sec_vec_.add( dof_indices_scalar_[0], normal_displacement_ );
        output_div_vec_.add( dof_indices_scalar_[0], normal_displacement_ / eq_fields_->cross_section(p_low) );
    }



private:
    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    FEValues<3> fv_;                                          ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fsv_;                                         ///< FEValues of side (neighbour integral) object

    LocDofVec dof_indices_;                                   ///< Vector of local DOF indices of vector fields
    LocDofVec dof_indices_scalar_;                            ///< Vector of local DOF indices of scalar fields
    LocDofVec dof_indices_tensor_;                            ///< Vector of local DOF indices of tensor fields
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in neighbour calculation.

    double normal_displacement_;                              ///< Holds constributions of normal displacement.
    arma::mat33 normal_stress_;                               ///< Holds constributions of normal stress.

    /// Data vectors of output fields (FieldFE).
    VectorMPI output_vec_;
    VectorMPI output_stress_vec_;
    VectorMPI output_von_mises_stress_vec_;
    VectorMPI output_mean_stress_vec_;
    VectorMPI output_cross_sec_vec_;
    VectorMPI output_div_vec_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};



/**
 * Container class for assembly of constraint matrix for contact condition.
 */
template <unsigned int dim>
class ConstraintAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "ConstraintAssemblyElasticity"; }

    /// Constructor.
    ConstraintAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::coupling;
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->cross_section_min;
    }

    /// Destructor.
    ~ConstraintAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_values_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_side_JxW_values | update_normal_vectors);

        n_dofs_ = fe_->n_dofs();
        dof_indices_.resize(n_dofs_);
        local_matrix_.resize(n_dofs_*n_dofs_);
        vec_view_side_ = &fe_values_side_.vector_view(0);
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        if (!cell_lower_dim.is_own()) return;
        
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
		dof_indices_ = cell_higher_dim.get_loc_dof_indices();

        vector<LongIdx> lower_nodes;
        for (unsigned int i=0; i<dim; i++)
            lower_nodes.push_back(cell_lower_dim.elm().node(i).idx());

        vector<double> ignore_dofs(n_dofs_, false);
        for (unsigned int idof=0; idof<n_dofs_; idof++) {
            if ( dof_indices_[idof] >= 0 && cell_higher_dim.cell_dof(idof).dim == 0 && 
                 find(lower_nodes.begin(), lower_nodes.end(), cell_higher_dim.elm().node(cell_higher_dim.cell_dof(idof).n_face_idx).idx()) == lower_nodes.end() )
            {
                ignore_dofs[idof] = true;
            }
        }

		fe_values_side_.reinit(neighb_side.side());

        for (unsigned int i=0; i<n_dofs_; i++)
            local_matrix_[i] = 0;

        // Assemble matrix and vector for contact conditions in the form B*x <= c,
        // where B*x is the average jump of normal displacements and c is the average cross-section on element.
        // Positive value means that the fracture closes.
        unsigned int k=0;
        double local_vector = 0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = fe_values_side_.normal_vector(k);

            local_vector += (eq_fields_->cross_section(p_low) - eq_fields_->cross_section_min(p_low))*fe_values_side_.JxW(k) / cell_lower_dim.elm().measure() / cell_lower_dim.elm()->n_neighs_vb();

            for (unsigned int i=0; i<n_dofs_; i++)
            {
                if (ignore_dofs[i]) continue;
                local_matrix_[i] += eq_fields_->cross_section(p_high)*arma::dot(vec_view_side_->value(i,k), nv)*fe_values_side_.JxW(k) / cell_lower_dim.elm().measure();
            }
        	k++;
        }

        int arow_local[1] = { eq_data_->constraint_idx_local[cell_lower_dim.elm_idx()] };
        MatSetValuesLocal(eq_data_->constraint_matrix, 1, arow_local, n_dofs_, dof_indices_.memptr(), &(local_matrix_[0]), ADD_VALUES);
        VecSetValueLocal(eq_data_->constraint_vec, arow_local[0], local_vector, ADD_VALUES);
    }



private:



    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    FEValues<3> fe_values_side_;                              ///< FEValues of side object

    LocDofVec dof_indices_;                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in boundary / neighbour calculation.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


void rref(arma::mat& M, double tol) {
    const unsigned int row_count = M.n_rows;
    const unsigned int column_count = M.n_cols;
    if (row_count == 0 || column_count <= 1) return;

    // The last column is the right-hand side; pivot only in coefficient columns.
    const unsigned int coefficient_columns = column_count - 1;
    unsigned int row = 0;

    for (unsigned int lead = 0; lead < coefficient_columns && row < row_count; ++lead) {
        unsigned int pivot_row = row;
        double pivot_abs = std::abs(M(pivot_row, lead));
        for (unsigned int i = row + 1; i < row_count; ++i) {
            const double value_abs = std::abs(M(i, lead));
            if (value_abs > pivot_abs) {
                pivot_abs = value_abs;
                pivot_row = i;
            }
        }

        if (pivot_abs <= tol) {
            M.submat(row, lead, row_count - 1, lead).zeros();
            continue;
        }

        M.swap_rows(row, pivot_row);
        M.row(row) /= M(row, lead);

        for (unsigned int i = 0; i < row_count; ++i) {
            if (i == row) continue;
            const double factor = M(i, lead);
            if (std::abs(factor) > tol)
                M.row(i) -= factor * M.row(row);
            else
                M(i, lead) = 0.0;
        }

        ++row;
    }

    M.elem(arma::find(arma::abs(M) <= tol)).zeros();
}


template <unsigned int dim>
class DirichletAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "DirichletAssemblyElasticity"; }

    /// Constructor.
    DirichletAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data), drop_tol(1e-12) {
        this->active_integrals_ = ActiveIntegrals::boundary;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_displacement;
    }

    /// Destructor.
    ~DirichletAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_values_bdr_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
        n_dofs_ = fe_->n_dofs();
        dof_indices_.resize(n_dofs_);
        vec_view_bdr_ = &fe_values_bdr_side_.vector_view(0);
    }


    /// Assembles boundary integral.
    void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dof_indices_ = dh_cell.get_loc_dof_indices();
        fe_values_bdr_side_.reinit(side);

        unsigned int bc_type;
        {
            auto p_side = *( this->boundary_points(cell_side).begin() );
    	    auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
    	    bc_type = eq_fields_->bc_type(p_bdr);
        }

        unsigned int k = 0;
        if (bc_type == EqFields::bc_type_displacement)
        {
            arma::mat A(n_dofs_,n_dofs_+1, arma::fill::zeros);
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                auto bc_displacement = eq_fields_->bc_displacement(p_bdr);
                for (unsigned int i=0; i<n_dofs_; i++) {
                    for (unsigned int j=0; j<n_dofs_; j++)
                        A(i,j) += arma::dot(vec_view_bdr_->value(j,k), vec_view_bdr_->value(i,k)) *
                            fe_values_bdr_side_.JxW(k);
                    A(i,n_dofs_) += arma::dot(bc_displacement, vec_view_bdr_->value(i,k)) *
                            fe_values_bdr_side_.JxW(k);
                }
                ++k;
            }
            rref(A, drop_tol);
            for (unsigned int row_idx=0; row_idx<A.n_rows; row_idx++) {
                arma::uvec indices = arma::find(arma::abs(A(row_idx, arma::span(0,A.n_cols-2))) > drop_tol);
                ASSERT_LE(indices.size(), 1).error("Dirichlet constraint applies to more than one dof!");
                if (indices.size() == 1) {
                    if (std::find(eq_data_->dirichlet_dofs.begin(), eq_data_->dirichlet_dofs.end(), dof_indices_[indices(0)]) == eq_data_->dirichlet_dofs.end()) {
                        eq_data_->dirichlet_dofs.push_back(dof_indices_[indices(0)]);
                        eq_data_->dirichlet_coefs.push_back(A(row_idx,indices(0)));
                        eq_data_->dirichlet_values.push_back(A(row_idx,A.n_cols-1));
                        eq_data_->dirichlet_row_starts.push_back( eq_data_->dirichlet_row_starts.back() + indices.size() );
                    }
                }
            }

        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            arma::mat A(n_dofs_,n_dofs_+1, arma::fill::zeros);
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                normal_vector_ = fe_values_bdr_side_.normal_vector(k);
                auto bc_displacement = eq_fields_->bc_displacement(p_bdr);
                for (unsigned int i=0; i<n_dofs_; i++) {
                    for (unsigned int j=0; j<n_dofs_; j++)
                        A(i,j) += arma::dot(vec_view_bdr_->value(j,k), normal_vector_) *
                            arma::dot(vec_view_bdr_->value(i,k), normal_vector_) *
                            fe_values_bdr_side_.JxW(k);
                    A(i,n_dofs_) += arma::dot(bc_displacement, normal_vector_) *
                            arma::dot(vec_view_bdr_->value(i,k),normal_vector_) *
                            fe_values_bdr_side_.JxW(k);
                }
                ++k;
            }
            rref(A, drop_tol);
            for (unsigned int row_idx=0; row_idx<A.n_rows; row_idx++) {
                arma::uvec indices = arma::find(arma::abs(A(row_idx, arma::span(0,A.n_cols-2))) > drop_tol);
                // ASSERT_LE(indices.size(), 1).error("Dirichlet constraint applies to more than one dof!");
                if (indices.size() > 0) {
                    for (unsigned int i=0; i<indices.size(); i++) {
                        eq_data_->dirichlet_dofs.push_back(dof_indices_[indices(i)]);
                        eq_data_->dirichlet_coefs.push_back(A(row_idx,indices(i)));
                    }
                    eq_data_->dirichlet_values.push_back(A(row_idx,A.n_cols-1));
                    eq_data_->dirichlet_row_starts.push_back( eq_data_->dirichlet_row_starts.back() + indices.size() );
                }
            }
        }

    }





private:
    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    FEValues<3> fe_values_bdr_side_;                          ///< FEValues of side (boundary integral) object

    LocDofVec dof_indices_;                             ///< Vector of global DOF indices
    const FEValuesViews::Vector<3> * vec_view_bdr_;           ///< Vector view in boundary calculation.

    arma::vec3 normal_vector_;

    const double drop_tol;                              // Tolerance to drop small values.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


class RigidBodyModesAssemblyElasticity
{
public:
    typedef Elasticity::EqData EqData;

    RigidBodyModesAssemblyElasticity(EqData *eq_data)
    : eq_data_(eq_data), max_dim_(0), n_modes_(0) {}

    ~RigidBodyModesAssemblyElasticity() {}

    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh)
    {
        max_dim_ = find_largest_element_dimension(dh);
        detect_components(dh);
        compute_component_modes(dh);
        print_component_diagnostics();
        create_modes_matrix(dh);

        for (auto cell : dh->own_range())
            assemble_cell(cell);
    }

private:
    unsigned int find_largest_element_dimension(std::shared_ptr<DOFHandlerMultiDim> dh) const
    {
        unsigned int max_dim = 0;
        for (auto cell : dh->local_range())
            max_dim = std::max(max_dim, cell.dim());
        return max_dim;
    }

    void detect_components(std::shared_ptr<DOFHandlerMultiDim> dh)
    {
        cell_component_.clear();
        component_origins_.clear();

        // Traverse the local cell graph in all dimensions. Components are seeded
        // only by max-dimensional cells, but lower-dimensional cells can connect
        // max-dimensional cells through fractures and mixed-dimensional joins.
        for (auto cell : dh->own_range())
            cell_component_[cell.elm_idx()] = -1;

        std::map<unsigned int, std::vector<unsigned int>> cell_neighbors;
        auto add_edge = [&](unsigned int a, unsigned int b) {
            if (a == b) return;
            if (cell_component_.find(a) == cell_component_.end()) return;
            if (cell_component_.find(b) == cell_component_.end()) return;

            cell_neighbors[a].push_back(b);
            cell_neighbors[b].push_back(a);
        };

        for (auto cell : dh->own_range()) {
            if (cell.dim() > 0) {
                for (DHCellSide side : cell.side_range()) {
                    for (DHCellSide edge_side : side.edge_sides())
                        add_edge(cell.elm_idx(), edge_side.cell().elm_idx());
                }
            }

            for (DHCellSide neighb_side : cell.neighb_sides())
                add_edge(cell.elm_idx(), neighb_side.cell().elm_idx());
        }

        int component_id = 0;
        for (auto cell : dh->own_range()) {
            if (cell.dim() != max_dim_) continue;
            if (cell_component_[cell.elm_idx()] >= 0) continue;

            std::queue<unsigned int> queue;
            queue.push(cell.elm_idx());
            cell_component_[cell.elm_idx()] = component_id;

            while (!queue.empty()) {
                unsigned int current_elm_idx = queue.front();
                queue.pop();

                for (unsigned int neighbor_elm_idx : cell_neighbors[current_elm_idx]) {
                    auto component_it = cell_component_.find(neighbor_elm_idx);
                    if (component_it == cell_component_.end()) continue;
                    if (component_it->second >= 0) continue;

                    component_it->second = component_id;
                    queue.push(neighbor_elm_idx);
                }
            }

            ++component_id;
        }
    }

    arma::vec3 normalized_vector(const arma::vec3 &v, const arma::vec3 &fallback) const
    {
        const double norm = arma::norm(v);
        if (norm > 0.0) return v / norm;
        return fallback;
    }

    arma::vec3 eigenvector_3d(const arma::mat &eigenvectors, unsigned int column) const
    {
        arma::vec3 v;
        v.zeros();
        for (unsigned int i=0; i<3; ++i)
            v[i] = eigenvectors(i, column);
        return v;
    }

    void init_component_normals(unsigned int component, const arma::mat33 &covariance)
    {
        component_normal_vectors_[component].clear();

        arma::vec eigenvalues;
        arma::mat eigenvectors;
        bool eig_ok = arma::eig_sym(eigenvalues, eigenvectors, covariance);

        arma::vec3 ex, ey, ez;
        ex.zeros();
        ey.zeros();
        ez.zeros();
        ex[0] = 1.0;
        ey[1] = 1.0;
        ez[2] = 1.0;

        if (max_dim_ == 3 || !eig_ok) {
            component_normal_vectors_[component].push_back(ex);
            component_normal_vectors_[component].push_back(ey);
            component_normal_vectors_[component].push_back(ez);
            return;
        }

        const double lambda_max = std::max(0.0, eigenvalues[eigenvalues.n_elem-1]);
        const double geometry_tol = 1e-8;

        if (max_dim_ == 2) {
            arma::vec3 normal = normalized_vector(eigenvector_3d(eigenvectors, 0), ez);

            if (lambda_max > 0.0) {
                const double relative_thickness =
                    std::sqrt(std::max(0.0, eigenvalues[0]) / lambda_max);
                if (relative_thickness > geometry_tol) {
                    PetscMPIInt rank;
                    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
                        "[feti user nullspace components] rank: %d, component: %d, "
                        "warning: 2D component is not planar, relative_thickness: %.16e\n",
                        rank, (int)component, relative_thickness);
                    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
                }
            }

            component_normal_vectors_[component].push_back(normal);
            return;
        }

        if (max_dim_ == 1) {
            arma::vec3 tangent = normalized_vector(eigenvector_3d(eigenvectors, 2), ex);

            if (lambda_max > 0.0) {
                const double relative_width =
                    std::sqrt(std::max(0.0, eigenvalues[1]) / lambda_max);
                if (relative_width > geometry_tol) {
                    PetscMPIInt rank;
                    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
                        "[feti user nullspace components] rank: %d, component: %d, "
                        "warning: 1D component is not linear, relative_width: %.16e\n",
                        rank, (int)component, relative_width);
                    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
                }
            }

            arma::vec3 helper = ex;
            if (std::abs(arma::dot(tangent, ey)) < std::abs(arma::dot(tangent, helper)))
                helper = ey;
            if (std::abs(arma::dot(tangent, ez)) < std::abs(arma::dot(tangent, helper)))
                helper = ez;

            arma::vec3 n1 = normalized_vector(arma::cross(tangent, helper), ey);
            arma::vec3 n2 = normalized_vector(arma::cross(tangent, n1), ez);
            component_normal_vectors_[component].push_back(n1);
            component_normal_vectors_[component].push_back(n2);
        }
    }

    void compute_component_modes(std::shared_ptr<DOFHandlerMultiDim> dh)
    {
        unsigned int n_components = 0;
        for (const auto &item : cell_component_) {
            if (item.second >= 0)
                n_components = std::max(n_components, static_cast<unsigned int>(item.second + 1));
        }

        component_origins_.assign(n_components, arma::vec3());
        component_normal_vectors_.assign(n_components, std::vector<arma::vec3>());
        component_first_mode_.assign(n_components, 0);
        component_n_modes_.assign(n_components, 0);

        std::vector<unsigned int> n_points(n_components, 0);
        for (auto &origin : component_origins_)
            origin.zeros();

        for (auto cell : dh->own_range()) {
            if (cell.dim() != max_dim_) continue;
            const int cell_component = component_of_cell(cell);
            if (cell_component < 0) continue;

            for (unsigned int i=0; i<cell.elm()->n_nodes(); ++i) {
                component_origins_[cell_component] += *cell.elm().node(i);
                ++n_points[cell_component];
            }
        }

        for (unsigned int c=0; c<component_origins_.size(); ++c) {
            if (n_points[c] > 0)
                component_origins_[c] /= n_points[c];
        }

        std::vector<arma::mat33> covariances(n_components);
        for (auto &covariance : covariances)
            covariance.zeros();

        for (auto cell : dh->own_range()) {
            if (cell.dim() != max_dim_) continue;
            const int cell_component = component_of_cell(cell);
            if (cell_component < 0) continue;

            for (unsigned int i=0; i<cell.elm()->n_nodes(); ++i) {
                arma::vec3 x = *cell.elm().node(i) - component_origins_[cell_component];
                covariances[cell_component] += x * x.t();
            }
        }

        n_modes_ = 0;
        for (unsigned int c=0; c<n_components; ++c) {
            init_component_normals(c, covariances[c]);
            component_first_mode_[c] = n_modes_;
            component_n_modes_[c] = 3 + component_normal_vectors_[c].size();
            n_modes_ += component_n_modes_[c];
        }
    }

    int component_of_cell(DHCellAccessor cell) const
    {
        auto component_it = cell_component_.find(cell.elm_idx());
        if (component_it != cell_component_.end() && component_it->second >= 0)
            return component_it->second;

        if (cell.dim() >= max_dim_) return -1;

        // Fallback for lower-dimensional local cells not reached during BFS.
        // Normally they are assigned by the all-dimensional traversal.
        int component = -1;
        for (DHCellSide neighb_side : cell.neighb_sides()) {
            DHCellAccessor neighbor = neighb_side.cell();
            if (neighbor.dim() != max_dim_) continue;

            auto neighbor_component_it = cell_component_.find(neighbor.elm_idx());
            if (neighbor_component_it == cell_component_.end()) continue;
            if (neighbor_component_it->second < 0) continue;

            if (component < 0) {
                component = neighbor_component_it->second;
            } else if (component != neighbor_component_it->second) {
                // Ambiguous lower-dimensional cell on an interface between local
                // components. Keep deterministic ownership to avoid writing one
                // DOF into multiple RBM blocks.
                component = std::min(component, neighbor_component_it->second);
            }
        }

        return component;
    }

    void create_modes_matrix(std::shared_ptr<DOFHandlerMultiDim> dh)
    {
        if (eq_data_->rigid_body_modes_matrix != nullptr)
            MatDestroy(&eq_data_->rigid_body_modes_matrix);

        PetscInt first_mode = 0;
        PetscInt n_modes_petsc = n_modes_;
        MPI_Exscan(&n_modes_petsc, &first_mode, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);

        std::vector<PetscInt> cols_l2g(n_modes_);
        for (unsigned int i=0; i<n_modes_; ++i)
            cols_l2g[i] = first_mode + i;

        ISLocalToGlobalMapping row_l2g = NULL;
        ISLocalToGlobalMapping col_l2g = NULL;

        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1,
            dh->get_local_to_global_map().size(),
            dh->get_local_to_global_map().data(),
            PETSC_COPY_VALUES, &row_l2g);

        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1,
            n_modes_, cols_l2g.data(), PETSC_COPY_VALUES, &col_l2g);

        MatCreateIS(PETSC_COMM_WORLD, 1,
                    dh->lsize(), n_modes_,
                    PETSC_DETERMINE, PETSC_DETERMINE,
                    row_l2g, col_l2g,
                    &eq_data_->rigid_body_modes_matrix);

        Mat local_modes = NULL;
        MatISGetLocalMat(eq_data_->rigid_body_modes_matrix, &local_modes);
        MatSeqAIJSetPreallocation(local_modes, n_modes_, NULL);
        MatISRestoreLocalMat(eq_data_->rigid_body_modes_matrix, &local_modes);

        ISLocalToGlobalMappingDestroy(&row_l2g);
        ISLocalToGlobalMappingDestroy(&col_l2g);

        MatZeroEntries(eq_data_->rigid_body_modes_matrix);
    }

    void print_component_diagnostics() const
    {
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "[feti user nullspace components] rank: %d, components: %d, "
            "local_modes: %d, component_modes:",
            rank, (int)component_origins_.size(), (int)n_modes_);
        for (unsigned int c=0; c<component_n_modes_.size(); ++c)
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, " %d", (int)component_n_modes_[c]);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }

    void assemble_cell(DHCellAccessor cell)
    {
        LocDofVec dof_indices = cell.get_loc_dof_indices();

        const int component = component_of_cell(cell);
        if (component < 0) return;

        const unsigned int first_mode = component_first_mode_[component];
        const arma::vec3 &origin = component_origins_[component];
        const std::vector<arma::vec3> &normal_vectors = component_normal_vectors_[component];

        for (unsigned int i=0; i<cell.n_dofs(); ++i) {
            const Dof &dof = cell.cell_dof(i);
            const arma::vec &dof_coefs = dof.coefs;
            arma::vec3 x = dof_real_point(cell, dof);
            arma::vec3 mode_value;

            for (unsigned int mode=0; mode<3; ++mode) {
                mode_value.zeros();
                mode_value[mode] = 1.0;
                MatSetValueLocal(eq_data_->rigid_body_modes_matrix,
                                 dof_indices[i], first_mode + mode, arma::dot(dof_coefs, mode_value), INSERT_VALUES);
            }

            for (unsigned int n=0; n<normal_vectors.size(); ++n) {
                mode_value = arma::cross(normal_vectors[n], x - origin);
                MatSetValueLocal(eq_data_->rigid_body_modes_matrix,
                                 dof_indices[i], first_mode + 3+n, arma::dot(dof_coefs, mode_value), INSERT_VALUES);
            }
        }
    }

    arma::vec3 dof_real_point(DHCellAccessor cell, const Dof &dof) const
    {
        arma::vec3 x;
        x.zeros();
        for (unsigned int i=0; i<dof.coords.n_elem; ++i)
            x += dof.coords[i] * (*cell.elm().node(i));
        return x;
    }

    EqData *eq_data_;

    unsigned int max_dim_;
    unsigned int n_modes_;

    std::vector<arma::vec3> component_origins_;
    std::vector<std::vector<arma::vec3>> component_normal_vectors_;
    std::vector<unsigned int> component_first_mode_;
    std::vector<unsigned int> component_n_modes_;
    std::map<unsigned int, int> cell_component_;
};



#endif /* ASSEMBLY_ELASTICITY_HH_ */
