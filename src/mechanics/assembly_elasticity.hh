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

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "mechanics/elasticity.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class StiffnessAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqData EqData;

    /// Constructor.
    StiffnessAssemblyElasticity(EqData *data)
    : AssemblyBase<dim>(1), data_(data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        std::vector<string> sub_names = {"X", "d", "lame_mu", "lame_lambda", "dirichlet_penalty", "young_modulus",
                            "poisson_ratio", "cross_section", "bc_type", "fracture_sigma"};
        this->used_fields_ = data_->subset( sub_names );
    }

    /// Destructor.
    ~StiffnessAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(FMT_UNUSED std::shared_ptr<Balance> balance) {
        //this->balance_ = balance;

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
    inline void cell_integral(unsigned int element_patch_idx, unsigned int dh_local_idx)
    {
        if ((int)dh_local_idx == -1) return;
        DHCellAccessor cell(data_->dh_.get(), dh_local_idx);
        if (cell.dim() != dim) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        unsigned int k=0;
        for (auto p : data_->stiffness_assembly_->bulk_points(element_patch_idx, cell.dim()) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += data_->cross_section(p)*(
                                                2*data_->lame_mu(p)*arma::dot(vec_view_->sym_grad(j,k), vec_view_->sym_grad(i,k))
                                                + data_->lame_lambda(p)*vec_view_->divergence(j,k)*vec_view_->divergence(i,k)
                                               )*fe_values_.JxW(k);
            }
            k++;
        }
        data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ_DBG(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);
        fe_values_side_.reinit(side);

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( data_->stiffness_assembly_->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = data_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqData::bc_type_displacement)
        {
            unsigned int k=0;
            for (auto p : data_->stiffness_assembly_->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (data_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k),vec_view_side_->value(j,k)) * fe_values_side_.JxW(k);
                k++;
            }
        }
        else if (bc_type == EqData::bc_type_displacement_normal)
        {
            unsigned int k=0;
            for (auto p : data_->stiffness_assembly_->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (data_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k), fe_values_side_.normal_vector(k)) *
                                arma::dot(vec_view_side_->value(j,k), fe_values_side_.normal_vector(k)) * fe_values_side_.JxW(k);
                k++;
            }
        }

        data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }


    /// Assembles between elements of different dimensions.
    inline void neigbour_integral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

		cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
		ElementAccessor<3> cell_sub = cell_lower_dim.elm();
		fe_values_sub_.reinit(cell_sub);

		DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
		cell_higher_dim.get_dof_indices(side_dof_indices_[1]);
		fe_values_side_.reinit(neighb_side.side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<n_dofs_; i++)
                for (unsigned int m=0; m<2; ++m)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_ngh_[n][m][i*(n_dofs_)+j] = 0;

        // set transmission conditions
        unsigned int k=0;
        for (auto p_high : data_->stiffness_assembly_->coupling_points(neighb_side) )
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
                    arma::mat33 gvft = (n==0) ? vec_view_sub_->grad(i,k) : arma::zeros(3,3);

                    for (int m=0; m<2; m++)
                        for (unsigned int j=0; j<n_dofs_ngh_[m]; j++) {
                            arma::vec3 ui = (m==0) ? arma::zeros(3) : vec_view_side_->value(j,k);
                            arma::vec3 uf = (m==1) ? arma::zeros(3) : vec_view_sub_->value(j,k);
                            arma::mat33 guit = (m==1) ? mat_t(vec_view_side_->grad(j,k),nv) : arma::zeros(3,3);
                            double divuit = (m==1) ? arma::trace(guit) : 0;

                            local_matrix_ngh_[n][m][i*n_dofs_ngh_[m] + j] +=
                                    data_->fracture_sigma(p_low)*(
                                     arma::dot(vf-vi,
                                      2/data_->cross_section(p_low)*(data_->lame_mu(p_low)*(uf-ui)+(data_->lame_mu(p_low)+data_->lame_lambda(p_low))*(arma::dot(uf-ui,nv)*nv))
                                      + data_->lame_mu(p_low)*arma::trans(guit)*nv
                                      + data_->lame_lambda(p_low)*divuit*nv
                                     )
                                     - arma::dot(gvft, data_->lame_mu(p_low)*arma::kron(nv,ui.t()) + data_->lame_lambda(p_low)*arma::dot(ui,nv)*arma::eye(3,3))
                                    )*fe_values_sub_.JxW(k);
                        }

                }
            }
        	k++;
        }

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int m=0; m<2; ++m)
                data_->ls->mat_set_values(n_dofs_ngh_[n], side_dof_indices_[n].data(), n_dofs_ngh_[m], side_dof_indices_[m].data(), &(local_matrix_ngh_[n][m][0]));
    }


    /// Implements @p AssemblyBase::reallocate_cache.
    void reallocate_cache(const ElementCacheMap &cache_map) override
    {
        used_fields_.set_dependency();
        used_fields_.cache_reallocate(cache_map);
    }


private:
    inline arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
    {
      arma::mat33 mt = m - m*arma::kron(n,n.t());
      return mt;
    }



    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data object shared with EqData
    EqData *data_;

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

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class RhsAssemblyElasticity : public AssemblyBase<dim>
{
public:
    typedef typename Elasticity::EqData EqData;

    /// Constructor.
    RhsAssemblyElasticity(EqData *data)
    : AssemblyBase<dim>(1), data_(data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        std::vector<string> sub_names = {"X", "d", "load", "potential_load", "cross_section", "dirichlet_penalty",
                            "lame_mu", "lame_lambda", "young_modulus", "poisson_ratio", "bc_type", "bc_displacement",
		                    "bc_traction", "fracture_sigma"};
        this->used_fields_ = data_->subset( sub_names );
    }

    /// Destructor.
    ~RhsAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(FMT_UNUSED std::shared_ptr<Balance> balance) {
        //this->balance_ = balance;

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
    inline void cell_integral(unsigned int element_patch_idx, unsigned int dh_local_idx)
    {
        if ((int)dh_local_idx == -1) return;
        DHCellAccessor cell(data_->dh_.get(), dh_local_idx);
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        //local_source_balance_vector.assign(n_dofs_, 0);
        //local_source_balance_rhs.assign(n_dofs_, 0);

        // compute sources
        unsigned int k=0;
        for (auto p : data_->rhs_assembly_->bulk_points(element_patch_idx, cell.dim()) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += (
                                 arma::dot(data_->load(p), vec_view_->value(i,k))
                                 -data_->potential_load(p)*vec_view_->divergence(i,k)
                                )*data_->cross_section(p)*fe_values_.JxW(k);
            ++k;
        }
        data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));

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
    	ASSERT_EQ_DBG(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);
        fe_values_bdr_side_.reinit(side);

        auto p_side = *( data_->rhs_assembly_->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = data_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        // local_flux_balance_vector.assign(n_dofs_, 0);
        // local_flux_balance_rhs = 0;

        unsigned int k=0;
        if (bc_type == EqData::bc_type_displacement)
        {
            double side_measure = cell_side.measure();
            for (auto p : data_->rhs_assembly_->boundary_points(cell_side) )
            {
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (data_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(data_->bc_displacement(p), vec_view_bdr_->value(i,k)) *
					        fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqData::bc_type_displacement_normal)
        {
            double side_measure = cell_side.measure();
            for (auto p : data_->rhs_assembly_->boundary_points(cell_side) )
            {
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (data_->dirichlet_penalty(p) / side_measure) *
                            arma::dot(data_->bc_displacement(p), fe_values_bdr_side_.normal_vector(k)) *
                            arma::dot(vec_view_bdr_->value(i,k), fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqData::bc_type_traction)
        {
            for (auto p : data_->rhs_assembly_->boundary_points(cell_side) )
            {
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += data_->cross_section(p) *
                            arma::dot(vec_view_bdr_->value(i,k), data_->bc_traction(p) + data_->potential_load(p) * fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));


//             balance_->add_flux_matrix_values(subst_idx, loc_b, side_dof_indices, local_flux_balance_vector);
//             balance_->add_flux_vec_value(subst_idx, loc_b, local_flux_balance_rhs);
		// ++loc_b;
    }


    /// Assembles between elements of different dimensions.
    inline void neigbour_integral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

		cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
		ElementAccessor<3> cell_sub = cell_lower_dim.elm();
		fe_values_sub_.reinit(cell_sub);

		DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
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
        for (auto p_high : data_->rhs_assembly_->coupling_points(neighb_side) )
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

                    local_rhs_ngh_[n][i] -= data_->fracture_sigma(p_low) * data_->cross_section(p_high) *
                            arma::dot(vf-vi, data_->potential_load(p_high) * nv) * fe_values_sub_.JxW(k);
                }
            }
            ++k;
        }

        for (unsigned int n=0; n<2; ++n)
            data_->ls->rhs_set_values(n_dofs_ngh_[n], side_dof_indices_[n].data(), &(local_rhs_ngh_[n][0]));
    }


    /// Implements @p AssemblyBase::reallocate_cache.
    void reallocate_cache(const ElementCacheMap &cache_map) override
    {
        used_fields_.set_dependency();
        used_fields_.cache_reallocate(cache_map);
    }


private:
    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data object shared with EqData
    EqData *data_;

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
    typedef typename Elasticity::EqData EqData;

    /// Constructor.
    OutpuFieldsAssemblyElasticity(EqData *data)
    : AssemblyBase<dim>(0), data_(data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling);
        std::vector<string> sub_names = {"X", "d", "cross_section", "lame_mu", "lame_lambda", "young_modulus", "poisson_ratio"};
        this->used_fields_ = data_->subset( sub_names );
    }

    /// Destructor.
    ~OutpuFieldsAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(FMT_UNUSED std::shared_ptr<Balance> balance) {
        //this->balance_ = balance;

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

        output_vec_ = data_->output_field_ptr->vec();
        output_stress_vec_ = data_->output_stress_ptr->vec();
        output_von_mises_stress_vec_ = data_->output_von_mises_stress_ptr->vec();
        output_cross_sec_vec_ = data_->output_cross_section_ptr->vec();
        output_div_vec_ = data_->output_div_ptr->vec();
    }


    /// Assemble integral over element
    inline void cell_integral(unsigned int element_patch_idx, unsigned int dh_local_idx)
    {
        if ((int)dh_local_idx == -1) return;
        DHCellAccessor cell(data_->dh_.get(), dh_local_idx);
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;
        DHCellAccessor cell_tensor = cell.cell_with_other_dh(data_->dh_tensor_.get());
        DHCellAccessor cell_scalar = cell.cell_with_other_dh(data_->dh_scalar_.get());

        auto elm = cell.elm();

        fv_.reinit(elm);
        dof_indices_        = cell.get_loc_dof_indices();
        dof_indices_scalar_ = cell_scalar.get_loc_dof_indices();
        dof_indices_tensor_ = cell_tensor.get_loc_dof_indices();

        auto p = *( data_->outout_fields_assembly_->bulk_points(element_patch_idx, cell.dim()).begin() );

        arma::mat33 stress = arma::zeros(3,3);
        double div = 0;
        for (unsigned int i=0; i<n_dofs_; i++)
        {
            stress += (2*data_->lame_mu(p)*vec_view_->sym_grad(i,0) + data_->lame_lambda(p)*vec_view_->divergence(i,0)*arma::eye(3,3))*output_vec_[dof_indices_[i]];
            div += vec_view_->divergence(i,0)*output_vec_[dof_indices_[i]];
        }

        arma::mat33 stress_dev = stress - arma::trace(stress)/3*arma::eye(3,3);
        double von_mises_stress = sqrt(1.5*arma::dot(stress_dev, stress_dev));
        output_div_vec_[dof_indices_scalar_[0]] += div;

        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<3; j++)
                output_stress_vec_[dof_indices_tensor_[i*3+j]] += stress(i,j);
        output_von_mises_stress_vec_[dof_indices_scalar_[0]] = von_mises_stress;

        output_cross_sec_vec_[dof_indices_scalar_[0]] += data_->cross_section(p);
    }


    /// Assembles between elements of different dimensions.
    inline void neigbour_integral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

//        cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
//        ElementAccessor<3> cell_sub = cell_lower_dim.elm();
//        fe_values_sub_.reinit(cell_sub);

//        DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
//        cell_higher_dim.get_dof_indices(side_dof_indices_[1]);
//        fe_values_side_.reinit(neighb_side.side());

    }


    /// Implements @p AssemblyBase::reallocate_cache.
    void reallocate_cache(const ElementCacheMap &cache_map) override
    {
        used_fields_.set_dependency();
        used_fields_.cache_reallocate(cache_map);
    }


private:
    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

    /// Data object shared with EqData
    EqData *data_;

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

    /// Data vectors of output fields (FieldFE).
    VectorMPI output_vec_;
    VectorMPI output_stress_vec_;
    VectorMPI output_von_mises_stress_vec_;
    VectorMPI output_cross_sec_vec_;
    VectorMPI output_div_vec_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ASSEMBLY_ELASTICITY_HH_ */

