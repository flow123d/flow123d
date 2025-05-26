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
#include "fem/op_factory.hh"
#include "fem/patch_op_impl.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class StiffnessAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssemblyElasticity"; }

    /// Constructor.
    StiffnessAssemblyElasticity(EqFields *eq_fields, EqData *eq_data, std::shared_ptr<EvalPoints> eval_points, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values, eval_points), eq_fields_(eq_fields), eq_data_(eq_data), // quad_order = 1
      JxW_( this->bulk_values().JxW() ),
      JxW_side_( this->side_values().JxW() ),
      normal_( this->side_values().normal_vector() ),
      deform_side_( this->side_values().vector_shape() ),
      grad_deform_( this->bulk_values().grad_vector_shape() ),
      sym_grad_deform_( this->bulk_values().vector_sym_grad() ),
      div_deform_( this->bulk_values().vector_divergence() ),
      deform_join_( this->join_values().vector_join_shape() ),
      deform_join_grad_( this->join_values().gradient_vector_join_shape() ) {
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

        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();
        n_dofs_sub_ = fe_low->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2*n_dofs_);
        local_matrix_.resize(4*n_dofs_*n_dofs_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        // Fracture stiffness is assembled in dimjoin_integral.
        if (cell.elm()->n_neighs_vb() > 0) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
                                                arma::dot(eq_fields_->stress_tensor(p,sym_grad_deform_.shape(j)(p)), sym_grad_deform_.shape(i)(p))
                                               )*JxW_(p);
            }
        }
        eq_data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqFields::bc_type_displacement)
        {
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p),deform_side_.shape(j)(p)) * JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p), normal_(p)) *
                                arma::dot(deform_side_.shape(j)(p), normal_(p)) * JxW_side_(p);
            }
        }

        eq_data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i] = dof_indices_[i];
        }

		DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i+n_dofs_ngh_[0]] = dof_indices_[i];
        }

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

		unsigned int n_neighs = cell_lower_dim.elm()->n_neighs_vb();

        for (unsigned int i=0; i<n_dofs_ngh_[0]+n_dofs_ngh_[1]; i++)
            for (unsigned int j=0; j<n_dofs_ngh_[0]+n_dofs_ngh_[1]; j++)
                local_matrix_[i*(n_dofs_ngh_[0]+n_dofs_ngh_[1])+j] = 0;

        // set transmission conditions
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;
                arma::vec3 diff_deform_i = deform_join_.shape(i)(p_low) - deform_join_.shape(i)(p_high);
                arma::mat33 grad_deform_i = deform_join_grad_.shape(i)(p_low);  // low dim element
                arma::mat33 semi_grad_i = grad_deform_i + n_neighs/eq_fields_->cross_section(p_low)*arma::kron(diff_deform_i,nv.t());
                arma::mat33 semi_sym_grad_i = 0.5*(semi_grad_i + semi_grad_i.t());

                for (uint j=0; j<deform_join_.n_dofs_both(); ++j) {
                    arma::vec3 deform_j_high = deform_join_.shape(j)(p_high);
                    arma::vec3 diff_deform_j = deform_join_.shape(j)(p_low) - deform_j_high;
                    arma::mat33 grad_deform_j = deform_join_grad_.shape(j)(p_low);  // low dim element
                    arma::mat33 semi_grad_j = grad_deform_j + n_neighs/eq_fields_->cross_section(p_low)*arma::kron(diff_deform_j,nv.t());
                    arma::mat33 semi_sym_grad_j = 0.5*(semi_grad_j + semi_grad_j.t());

                    local_matrix_[i * (n_dofs_ngh_[0]+n_dofs_ngh_[1]) + j] +=
                            (
                                     eq_fields_->fracture_sigma(p_low)*eq_fields_->cross_section(p_low) / n_neighs
                                     * arma::dot(semi_sym_grad_i, eq_fields_->stress_tensor(p_low,semi_sym_grad_j))

                                     // This term corrects the tangential part of the above product so that there is no
                                     // dependence on fracture_sigma.
                                     // TODO: Fracture_sigma should be possibly removed and replaced by anisotropic elasticity.
                                     + (1-eq_fields_->fracture_sigma(p_low))*eq_fields_->cross_section(p_low) / n_neighs
                                       * arma::dot(0.5*(grad_deform_i+grad_deform_i.t()), eq_fields_->stress_tensor(p_low,0.5*(grad_deform_j+grad_deform_j.t())))
                            )*JxW_side_(p_high);
                }
            }

        }

        eq_data_->ls->mat_set_values(n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), &(local_matrix_[0]));
    }



private:
    inline arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
    {
      arma::mat33 mt = m - m*arma::kron(n,n.t());
      return mt;
    }


    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                               ///< Number of dofs
    unsigned int n_dofs_sub_;                                           ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                              ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    vector<LongIdx> dof_indices_;                                       ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                                  ///< vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                                  ///< Auxiliary vector for assemble methods

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Tensor> sym_grad_deform_;
    FeQArray<Scalar> div_deform_;
    FeQJoin<Vector> deform_join_;
    FeQJoin<Tensor> deform_join_grad_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class RhsAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "RhsAssemblyElasticity"; }

    /// Constructor.
    RhsAssemblyElasticity(EqFields *eq_fields, EqData *eq_data, std::shared_ptr<EvalPoints> eval_points, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values, eval_points), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_( this->bulk_values().JxW() ),
      JxW_side_( this->side_values().JxW() ),
      normal_( this->side_values().normal_vector() ),
      deform_( this->bulk_values().vector_shape() ),
      deform_side_( this->side_values().vector_shape() ),
	  grad_deform_( this->bulk_values().grad_vector_shape() ),
      div_deform_( this->bulk_values().vector_divergence() ),
      deform_join_( this->join_values().vector_join_shape() ) {
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

        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();
        n_dofs_sub_ = fe_low->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(n_dofs_sub_ + n_dofs_);
        local_rhs_.resize(2*n_dofs_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        //local_source_balance_vector.assign(n_dofs_, 0);
        //local_source_balance_rhs.assign(n_dofs_, 0);

        // compute sources
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += (
                                 arma::dot(eq_fields_->load(p), deform_.shape(i)(p))
                                 -eq_fields_->potential_load(p)*div_deform_.shape(i)(p)
                                 -arma::dot(eq_fields_->initial_stress(p), grad_deform_.shape(i)(p))
                                )*eq_fields_->cross_section(p)*JxW_(p);
        }
        eq_data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));

//         for (unsigned int i=0; i<n_dofs_; i++)
//         {
//             for (unsigned int k=0; k<qsize_; k++) // point range
//                 local_source_balance_vector[i] -= 0;//sources_sigma[k]*fe_vals_[vec_view_].value(i,k)*JxW_(k);
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
        if (!cell_side.cell().is_own()) return;

        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        // local_flux_balance_vector.assign(n_dofs_, 0);
        // local_flux_balance_rhs = 0;

        // addtion from initial stress
        for (auto p : this->boundary_points(cell_side) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += eq_fields_->cross_section(p) *
                        arma::dot(( eq_fields_->initial_stress(p) * normal_(p)),
                                    deform_side_.shape(i)(p)) *
                        JxW_side_(p);
        }

        if (bc_type == EqFields::bc_type_displacement)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(eq_fields_->bc_displacement(p_bdr), deform_side_.shape(i)(p)) *
					        JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                            arma::dot(eq_fields_->bc_displacement(p_bdr), normal_(p)) *
                            arma::dot(deform_side_.shape(i)(p), normal_(p)) *
                            JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_traction)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(deform_side_.shape(i)(p), eq_fields_->bc_traction(p_bdr) + eq_fields_->ref_potential_load(p) * normal_(p)) *
                            JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_stress)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    // stress is multiplied by inward normal to obtain traction
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(deform_side_.shape(i)(p), -eq_fields_->bc_stress(p_bdr)*normal_(p)
                            + eq_fields_->ref_potential_load(p) * normal_(p))
                            * JxW_side_(p);
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

        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i] = dof_indices_[i];
        }

	    DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i+n_dofs_ngh_[0]] = dof_indices_[i];
        }

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int i=0; i<2*n_dofs_; i++)
            local_rhs_[i] = 0;

        // set transmission conditions
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;

                arma::vec3 vi = deform_join_.shape(i)(p_high);
                arma::vec3 vf = deform_join_.shape(i)(p_low);

                local_rhs_[i] -= eq_fields_->fracture_sigma(p_low) * eq_fields_->cross_section(p_high) *
                        arma::dot(vf-vi, eq_fields_->potential_load(p_high) * nv) * JxW_side_(p_high);
            }
        }

        eq_data_->ls->rhs_set_values(n_dofs_ngh_[0]+n_dofs_ngh_[1], side_dof_indices_.data(), &(local_rhs_[0]));
    }



private:
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                               ///< Number of dofs
    unsigned int n_dofs_sub_;                                           ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                              ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    vector<LongIdx> dof_indices_;                                       ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                                  ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_rhs_;                                     ///< Auxiliary vector for assemble methods

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Scalar> div_deform_;
    FeQJoin<Vector> deform_join_;


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

template <unsigned int dim>
class OutpuFieldsAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::OutputEqData EqData;

    static constexpr const char * name() { return "OutpuFieldsAssemblyElasticity"; }

    /// Constructor.
    OutpuFieldsAssemblyElasticity(EqFields *eq_fields, EqData *eq_data, std::shared_ptr<EvalPoints> eval_points, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values, eval_points), eq_fields_(eq_fields), eq_data_(eq_data),
      normal_( this->side_values().normal_vector() ),
      deform_side_( this->side_values().vector_shape() ),
	  grad_deform_( this->bulk_values().grad_vector_shape() ),
      sym_grad_deform_( this->bulk_values().vector_sym_grad() ),
      div_deform_( this->bulk_values().vector_divergence() ) {
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

        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();

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

        dof_indices_        = cell.get_loc_dof_indices();
        dof_indices_scalar_ = cell_scalar.get_loc_dof_indices();
        dof_indices_tensor_ = cell_tensor.get_loc_dof_indices();

        auto p = *( this->bulk_points(element_patch_idx).begin() );

        arma::mat33 stress = eq_fields_->initial_stress(p);
        double div = 0;
        for (unsigned int i=0; i<n_dofs_; i++)
        {
            stress += eq_fields_->stress_tensor(p,sym_grad_deform_.shape(i)(p))
                    * output_vec_.get(dof_indices_[i]);
            div += div_deform_.shape(i)(p)*output_vec_.get(dof_indices_[i]);
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

        dof_indices_ = cell_higher_dim.get_loc_dof_indices();
        auto p_high = *( this->coupling_points(neighb_side).begin() );
        auto p_low = p_high.lower_dim(cell_lower_dim);

        for (unsigned int i=0; i<n_dofs_; i++)
        {
            normal_displacement_ -= arma::dot(deform_side_.shape(i)(p_high)*output_vec_.get(dof_indices_[i]), normal_(p_high));
            arma::mat33 grad = -arma::kron(deform_side_.shape(i)(p_high)*output_vec_.get(dof_indices_[i]), normal_(p_high).t()) / eq_fields_->cross_section(p_low);
            normal_stress_ += eq_fields_->stress_tensor(p_low, 0.5*(grad+grad.t()));
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
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                               ///< Number of dofs
    LocDofVec dof_indices_;                                             ///< Vector of local DOF indices of vector fields
    LocDofVec dof_indices_scalar_;                                      ///< Vector of local DOF indices of scalar fields
    LocDofVec dof_indices_tensor_;                                      ///< Vector of local DOF indices of tensor fields

    double normal_displacement_;                                        ///< Holds constributions of normal displacement.
    arma::mat33 normal_stress_;                                         ///< Holds constributions of normal stress.

    /// Following data members represent Element quantities and FE quantities
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Tensor> sym_grad_deform_;
    FeQArray<Scalar> div_deform_;

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
class ConstraintAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename Elasticity::EqFields EqFields;
    typedef typename Elasticity::EqData EqData;

    static constexpr const char * name() { return "ConstraintAssemblyElasticity"; }

    /// Constructor.
    ConstraintAssemblyElasticity(EqFields *eq_fields, EqData *eq_data, std::shared_ptr<EvalPoints> eval_points, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values, eval_points), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_side_( this->side_values().JxW() ),
      normal_( this->side_values().normal_vector() ),
      deform_side_( this->side_values().vector_shape() ) {
        this->active_integrals_ = ActiveIntegrals::coupling;
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->cross_section_min;
    }

    /// Destructor.
    ~ConstraintAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();
        dof_indices_.resize(n_dofs_);
        local_matrix_.resize(n_dofs_*n_dofs_);
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        if (!cell_lower_dim.is_own()) return;
        
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        cell_higher_dim.get_dof_indices(dof_indices_);

        for (unsigned int i=0; i<n_dofs_; i++)
            local_matrix_[i] = 0;

        // Assemble matrix and vector for contact conditions in the form B*x <= c,
        // where B*x is the average jump of normal displacements and c is the average cross-section on element.
        // Positive value means that the fracture closes.
        double local_vector = 0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_(p_high);

            local_vector += (eq_fields_->cross_section(p_low) - eq_fields_->cross_section_min(p_low))*JxW_side_(p_high) / cell_lower_dim.elm().measure() / cell_lower_dim.elm()->n_neighs_vb();

            for (unsigned int i=0; i<n_dofs_; i++)
            {
                local_matrix_[i] += eq_fields_->cross_section(p_high)*arma::dot(deform_side_.shape(i)(p_high), nv)*JxW_side_(p_high) / cell_lower_dim.elm().measure();
            }
        }

        int arow[1] = { eq_data_->constraint_idx[cell_lower_dim.elm_idx()] };
        MatSetValues(eq_data_->constraint_matrix, 1, arow, n_dofs_, dof_indices_.data(), &(local_matrix_[0]), ADD_VALUES);
        VecSetValue(eq_data_->constraint_vec, arow[0], local_vector, ADD_VALUES);
    }



private:


    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                               ///< Number of dofs
    vector<LongIdx> dof_indices_;                                       ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;                         ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                                  ///< Auxiliary vector for assemble methods

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_side_;


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ASSEMBLY_ELASTICITY_HH_ */

