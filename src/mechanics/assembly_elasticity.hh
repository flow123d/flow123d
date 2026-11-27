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
#include "fem/patch_op_impl.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"
#include "la/linsys_PERMON.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class TEqData>
class StiffnessAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "Elasticity_Stiffness_Assembly"; }

    /// Constructor.
    StiffnessAssemblyElasticity(EqData *eq_data, PatchInternals *patch_internals)
    : AssemblyBasePatch<dim>(eq_data->quad_order(), patch_internals), eq_fields_(eq_data->eq_fields_.get()), eq_data_(eq_data), // quad_order = 1
      bulk_integral_( this->create_bulk_integral(this->quad_)),
      bdr_integral_( this->create_boundary_integral(this->quad_low_) ),
      coupling_integral_( this->create_coupling_integral(this->quad_) ),
      JxW_( bulk_integral_->JxW() ),
      JxW_side_( bdr_integral_->JxW() ),
      JxW_join_( coupling_integral_->JxW() ),
      normal_( bdr_integral_->normal_vector() ),
      normal_join_( coupling_integral_->normal_vector() ),
      deform_side_( bdr_integral_->vector_shape() ),
      grad_deform_( bulk_integral_->grad_vector_shape() ),
      sym_grad_deform_( bulk_integral_->vector_sym_grad() ),
      deform_join_( coupling_integral_->vector_join_shape() ),
      deform_join_grad_( coupling_integral_->gradient_vector_join_shape() ) {
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
    void initialize() {
        //this->balance_ = eq_data_->balance_;
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);

        n_dofs_ = this->n_dofs();
        n_dofs_high_ = this->n_dofs_high();
        n_dofs_ngh_ = { n_dofs_, n_dofs_high_ };
        dof_indices_.resize(n_dofs_high_);
        side_dof_indices_.resize(2*n_dofs_high_);
        local_matrix_.resize(4*n_dofs_high_*n_dofs_high_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        // Fracture stiffness is assembled in dimjoin_integral.
        if (cell.elm()->n_neighs_vb() > 0) return;

        dof_indices_ = cell.get_loc_dof_indices();

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        for (auto p : bulk_integral_->points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
                                                arma::dot(eq_fields_->stress_tensor(p,sym_grad_deform_.shape(j)(p)), sym_grad_deform_.shape(i)(p))
                                               )*JxW_(p);
            }
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

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( bdr_integral_->points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqFields::bc_type_displacement)
        {
            for (auto p : bdr_integral_->points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p),deform_side_.shape(j)(p)) * JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            for (auto p : bdr_integral_->points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p), normal_(p)) *
                                arma::dot(deform_side_.shape(j)(p), normal_(p)) * JxW_side_(p);
            }
        }

        eq_data_->ls->mat_set_values_local(n_dofs_, dof_indices_.memptr(), n_dofs_, dof_indices_.memptr(), &(local_matrix_[0]));
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 3) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim).error("Dimension of element mismatch!");

        dof_indices_ = cell_lower_dim.get_loc_dof_indices();
        for(unsigned int i=0; i<dof_indices_.n_elem; ++i) {
            side_dof_indices_[i] = dof_indices_[i];
        }

		DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        dof_indices_ = cell_higher_dim.get_loc_dof_indices();
        for(unsigned int i=0; i<dof_indices_.n_elem; ++i) {
            side_dof_indices_[i+n_dofs_ngh_[0]] = dof_indices_[i];
        }

        // TODO: Remove this unused code.
        // // mark dofs that do not contribute to dimjoin coupling
        // vector<LongIdx> lower_nodes;
        // for (unsigned int i=0; i<dim; i++)
        //     lower_nodes.push_back(cell_lower_dim.elm().node(i).idx());
        // vector<double> ignore_dofs(n_dofs_, false);

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

		unsigned int n_neighs = cell_lower_dim.elm()->n_neighs_vb();

        for (unsigned int i=0; i<n_dofs_ngh_[0]+n_dofs_ngh_[1]; i++)
            for (unsigned int j=0; j<n_dofs_ngh_[0]+n_dofs_ngh_[1]; j++)
                local_matrix_[i*(n_dofs_ngh_[0]+n_dofs_ngh_[1])+j] = 0;

        // set transmission conditions
        for (auto p_high : coupling_integral_->points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_join_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;
                arma::vec3 deform_i_low = deform_join_.shape(i)(p_low);
                arma::vec3 deform_i_high = deform_join_.shape(i)(p_high);
                arma::mat33 grad_deform_i = mat_t(deform_join_grad_.shape(i)(p_low),nv);  // low dim element
                if (eq_data_->fix_nullspace) {
                    arma::mat33 grad_deform_high = deform_join_grad_.shape(i)(p_high);
                    deform_i_high = deform_i_high - eq_fields_->cross_section(p_low)/2 * (grad_deform_high*nv);
                }
                arma::vec3 diff_deform_i = deform_i_low - deform_i_high;
                // TODO: Add the following line together with the commented out block of code below, implementing local_matrix_ngh assembly.
                // double div_grad_deform_i = arma::trace(grad_deform_i);
                arma::mat33 semi_grad_i = grad_deform_i + n_neighs/eq_fields_->cross_section(p_low)*arma::kron(diff_deform_i,nv.t());
                arma::mat33 semi_sym_grad_i = 0.5*(semi_grad_i + semi_grad_i.t());

                for (uint j=0; j<deform_join_.n_dofs_both(); ++j) {
                    arma::vec3 deform_j_low = deform_join_.shape(j)(p_low);
                    arma::vec3 deform_j_high = deform_join_.shape(j)(p_high);
                    arma::mat33 grad_deform_j = deform_join_grad_.shape(j)(p_low);  // low dim element
                    if (eq_data_->fix_nullspace) {
                        arma::mat33 grad_deform_high = deform_join_grad_.shape(j)(p_high);
                        deform_j_high = deform_j_high - eq_fields_->cross_section(p_low)/2 * (grad_deform_high*nv);
                    }
                    arma::vec3 diff_deform_j = deform_j_low - deform_j_high;
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
                            )*JxW_join_(p_high);
                    // TODO:
                    // - in old code:
                    // eq_fields_->fracture_sigma(p_low)*(
                    //  arma::dot(vf-vi,
                    //   2/eq_fields_->cross_section(p_low)*(eq_fields_->lame_mu(p_low)*(uf-ui)+(eq_fields_->lame_mu(p_low)+eq_fields_->lame_lambda(p_low))*(arma::dot(uf-ui,nv)*nv))
                    //   + eq_fields_->lame_mu(p_low)*arma::trans(guft)*nv
                    //   + eq_fields_->lame_lambda(p_low)*divuft*nv
                    //  )
                    //  - arma::dot(gvft, eq_fields_->lame_mu(p_low)*arma::kron(nv,ui.t()) + eq_fields_->lame_lambda(p_low)*arma::dot(ui,nv)*arma::eye(3,3))
                    // )*fe_values_sub_.JxW(k);
                    // - replace the above 7 lines starting "arma::dot(vf-vi," with the following code.
                    // It is mathematically equivalent, but the below code is more readable.
                    //     2/eq_fields_->cross_section(p_low)*(
                    //         eq_fields_->lame_mu(p_low)*arma::dot(vf-vi,uf-ui)
                    //       +(eq_fields_->lame_mu(p_low)+eq_fields_->lame_lambda(p_low))*arma::dot(uf-ui,nv)*arma::dot(vf-vi,nv)
                    //     )
                    //     + eq_fields_->lame_mu(p_low)*( arma::dot(vf-vi,guft.t()*nv) + arma::dot(uf-ui,gvft.t()*nv) )
                    //     + eq_fields_->lame_lambda(p_low)*( divuft*arma::dot(vf-vi,nv) + divvft*arma::dot(uf-ui,nv) ) // divvft > div_grad_deform_i
                    // )*fe_values_sub_.JxW(k) * 2/cell_lower_dim.elm()->n_neighs_vb();
                    // - use FeQJoin implementation
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
    unsigned int n_dofs_high_;                                          ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                              ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    LocDofVec dof_indices_;                                             ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                                  ///< vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                                  ///< Auxiliary vector for assemble methods

    std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;               ///< Bulk integral of assembly class
    std::shared_ptr<BoundaryIntegralAcc<dim>> bdr_integral_;            ///< Boundary integral of assembly class
    std::shared_ptr<CouplingIntegralAcc<dim>> coupling_integral_;       ///< Coupling integral of assembly class

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    FeQ<Scalar> JxW_join_;
    ElQ<Vector> normal_;
    ElQ<Vector> normal_join_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Tensor> sym_grad_deform_;
    FeQJoin<Vector> deform_join_;
    FeQJoin<Tensor> deform_join_grad_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim, class TEqData>
class RhsAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "Elasticity_Rhs_Assembly"; }

    /// Constructor.
    RhsAssemblyElasticity(EqData *eq_data, PatchInternals *patch_internals)
    : AssemblyBasePatch<dim>(eq_data->quad_order(), patch_internals), eq_fields_(eq_data->eq_fields_.get()), eq_data_(eq_data),
      bulk_integral_( this->create_bulk_integral(this->quad_)),
      bdr_integral_( this->create_boundary_integral(this->quad_low_) ),
      coupling_integral_( this->create_coupling_integral(this->quad_) ),
      JxW_( bulk_integral_->JxW() ),
      JxW_side_( bdr_integral_->JxW() ),
      JxW_join_( coupling_integral_->JxW() ),
      normal_( bdr_integral_->normal_vector() ),
      normal_join_( coupling_integral_->normal_vector() ),
      deform_( bulk_integral_->vector_shape() ),
      deform_side_( bdr_integral_->vector_shape() ),
	  grad_deform_( bulk_integral_->grad_vector_shape() ),
      div_deform_( bulk_integral_->vector_divergence() ),
      deform_join_( coupling_integral_->vector_join_shape() ) {
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
    void initialize() {
        //this->balance_ = eq_data_->balance_;
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);

        n_dofs_ = this->n_dofs();
        n_dofs_high_ = this->n_dofs_high();
        n_dofs_ngh_ = { n_dofs_, n_dofs_high_ };
        dof_indices_.resize(n_dofs_high_);
        side_dof_indices_.resize(n_dofs_ + n_dofs_high_);
        local_rhs_.resize(2*n_dofs_high_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        //local_source_balance_vector.assign(n_dofs_, 0);
        //local_source_balance_rhs.assign(n_dofs_, 0);

        // compute sources
        for (auto p : bulk_integral_->points(element_patch_idx) )
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

        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);

        auto p_side = *( bdr_integral_->points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        // local_flux_balance_vector.assign(n_dofs_, 0);
        // local_flux_balance_rhs = 0;

        // addtion from initial stress
        for (auto p : bdr_integral_->points(cell_side) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += eq_fields_->cross_section(p) *
                        arma::dot(( eq_fields_->initial_stress(p) * normal_(p)),
                                    deform_side_.shape(i)(p)) *
                        JxW_side_(p);
        }

        if (bc_type == EqFields::bc_type_displacement && !eq_data_->dirichlet_by_eq)
        {
            double side_measure = cell_side.measure();
            for (auto p : bdr_integral_->points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(eq_fields_->bc_displacement(p_bdr), deform_side_.shape(i)(p)) *
					        JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal && !eq_data_->dirichlet_by_eq)
        {
            double side_measure = cell_side.measure();
            for (auto p : bdr_integral_->points(cell_side) )
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
            for (auto p : bdr_integral_->points(cell_side) )
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
            for (auto p : bdr_integral_->points(cell_side) )
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
    	if (dim == 3) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim).error("Dimension of element mismatch!");

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
        for (auto p_high : coupling_integral_->points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_join_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;

                arma::vec3 vi = deform_join_.shape(i)(p_high);
                arma::vec3 vf = deform_join_.shape(i)(p_low);

                local_rhs_[i] -= eq_fields_->fracture_sigma(p_low) * eq_fields_->cross_section(p_high) *
                        arma::dot(vf-vi, eq_fields_->potential_load(p_high) * nv) * JxW_join_(p_high);
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
    unsigned int n_dofs_high_;                                          ///< Number of dofs (on higher dim element)
    std::vector<unsigned int> n_dofs_ngh_;                              ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    vector<LongIdx> dof_indices_;                                       ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                                  ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_rhs_;                                     ///< Auxiliary vector for assemble methods

    std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;               ///< Bulk integral of assembly class
    std::shared_ptr<BoundaryIntegralAcc<dim>> bdr_integral_;            ///< Boundary integral of assembly class
    std::shared_ptr<CouplingIntegralAcc<dim>> coupling_integral_;       ///< Coupling integral of assembly class

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    FeQ<Scalar> JxW_join_;
    ElQ<Vector> normal_;
    ElQ<Vector> normal_join_;
    FeQArray<Vector> deform_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Scalar> div_deform_;
    FeQJoin<Vector> deform_join_;


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

template <unsigned int dim, class TEqData>
class OutpuFieldsAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "Elasticity_OutpuFields_Assembly"; }

    /// Constructor.
    OutpuFieldsAssemblyElasticity(EqData *eq_data, PatchInternals *patch_internals)
    : AssemblyBasePatch<dim>(eq_data->quad_order(), patch_internals), eq_fields_(eq_data->eq_fields_.get()), eq_data_(eq_data),
      bulk_integral_( this->create_bulk_integral(this->quad_)),
      coupling_integral_( this->create_coupling_integral(this->quad_) ),
      normal_join_( coupling_integral_->normal_vector() ),
      deform_join_( coupling_integral_->vector_shape() ),
      grad_deform_( bulk_integral_->grad_vector_shape() ),
      sym_grad_deform_( bulk_integral_->vector_sym_grad() ),
      div_deform_( bulk_integral_->vector_divergence() ) {
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->lame_mu;
        this->used_fields_ += eq_fields_->lame_lambda;
        this->used_fields_ += eq_fields_->initial_stress;
    }

    /// Destructor.
    ~OutpuFieldsAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize() {
        //this->balance_ = eq_data_->balance_;

        n_dofs_ = this->n_dofs();
        n_dofs_high_ = this->n_dofs_high();

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

        auto p = *( bulk_integral_->points(element_patch_idx).begin() );

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
        if (dim == 3) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim).error("Dimension of element mismatch!");

        normal_displacement_ = 0;
        normal_stress_.zeros();

        DHCellAccessor cell_higher_dim = neighb_side.cell();
        DHCellAccessor cell_tensor = cell_lower_dim.cell_with_other_dh(eq_data_->dh_tensor_.get());
        DHCellAccessor cell_scalar = cell_lower_dim.cell_with_other_dh(eq_data_->dh_scalar_.get());

        dof_indices_ = cell_higher_dim.get_loc_dof_indices();
        auto p_high = *( coupling_integral_->points(neighb_side).begin() );
        auto p_low = p_high.lower_dim(cell_lower_dim);

        for (unsigned int i=0; i<n_dofs_high_; i++)
        {
            normal_displacement_ -= arma::dot(deform_join_.shape(i)(p_high)*output_vec_.get(dof_indices_[i]), normal_join_(p_high));
            arma::mat33 grad = -arma::kron(deform_join_.shape(i)(p_high)*output_vec_.get(dof_indices_[i]), normal_join_(p_high).t()) / eq_fields_->cross_section(p_low);
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
    unsigned int n_dofs_high_;                                          ///< Number of dofs of higher dim element
    LocDofVec dof_indices_;                                             ///< Vector of local DOF indices of vector fields
    LocDofVec dof_indices_scalar_;                                      ///< Vector of local DOF indices of scalar fields
    LocDofVec dof_indices_tensor_;                                      ///< Vector of local DOF indices of tensor fields

    double normal_displacement_;                                        ///< Holds constributions of normal displacement.
    arma::mat33 normal_stress_;                                         ///< Holds constributions of normal stress.

    std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;               ///< Bulk integral of assembly class
    std::shared_ptr<CouplingIntegralAcc<dim>> coupling_integral_;       ///< Coupling integral of assembly class

    /// Following data members represent Element quantities and FE quantities
    ElQ<Vector> normal_join_;
    FeQArray<Vector> deform_join_;
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
template <unsigned int dim, class TEqData>
class ConstraintAssemblyElasticity : public AssemblyBasePatch<dim>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "Elasticity_Constraint_Assembly"; }

    /// Constructor.
    ConstraintAssemblyElasticity(EqData *eq_data, PatchInternals *patch_internals)
    : AssemblyBasePatch<dim>(eq_data->quad_order(), patch_internals), eq_fields_(eq_data->eq_fields_.get()), eq_data_(eq_data),
      coupling_integral_( this->create_coupling_integral(this->quad_) ),
      JxW_join_( coupling_integral_->JxW() ),
      normal_join_( coupling_integral_->normal_vector() ),
      deform_join_( coupling_integral_->vector_shape() ) {
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->cross_section_min;
    }

    /// Destructor.
    ~ConstraintAssemblyElasticity() {}

    /// Initialize auxiliary vectors and other data members
    void initialize() {
        n_dofs_ = this->n_dofs_high();
        dof_indices_.resize(n_dofs_);
        local_matrix_.resize(n_dofs_*n_dofs_);
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 3) return;
        if (!cell_lower_dim.is_own()) return;
        
        ASSERT_EQ(cell_lower_dim.dim(), dim).error("Dimension of element mismatch!");

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

        for (unsigned int i=0; i<n_dofs_; i++)
            local_matrix_[i] = 0;

        // Assemble matrix and vector for contact conditions in the form B*x <= c,
        // where B*x is the average jump of normal displacements and c is the average cross-section on element.
        // Positive value means that the fracture closes.
        double local_vector = 0;
        for (auto p_high : coupling_integral_->points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_join_(p_high);

            local_vector += (eq_fields_->cross_section(p_low) - eq_fields_->cross_section_min(p_low))*JxW_join_(p_high) / cell_lower_dim.elm().measure() / cell_lower_dim.elm()->n_neighs_vb();

            for (unsigned int i=0; i<n_dofs_; i++)
            {
                if (ignore_dofs[i]) continue;
                local_matrix_[i] += eq_fields_->cross_section(p_high)*arma::dot(deform_join_.shape(i)(p_high), nv)*JxW_join_(p_high) / cell_lower_dim.elm().measure();
            }
        }

        int arow_local[1] = { eq_data_->constraint_idx_local[cell_lower_dim.elm_idx()] };
        MatSetValuesLocal(eq_data_->constraint_matrix, 1, arow_local, n_dofs_, dof_indices_.memptr(), &(local_matrix_[0]), ADD_VALUES);
        VecSetValueLocal(eq_data_->constraint_vec, arow_local[0], local_vector, ADD_VALUES);
    }



private:


    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                               ///< Number of dofs
    LocDofVec dof_indices_;                                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;                         ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                                  ///< Auxiliary vector for assemble methods

    std::shared_ptr<CouplingIntegralAcc<dim>> coupling_integral_;       ///< Coupling integral of assembly class

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_join_;
    ElQ<Vector> normal_join_;
    FeQArray<Vector> deform_join_;


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


void rref(arma::mat& M, double tol) {
    int lead = 0;
    int rowCount = M.n_rows;
    int columnCount = M.n_cols;

    for (int r = 0; r < rowCount; ++r) {
        if (columnCount <= lead) return;
        int i = r;
        while (fabs(M(i, lead)) < tol) {
            i++;
            if (rowCount == i) {
                i = r;
                lead++;
                if (columnCount == lead) return;
            }
        }
        M.swap_rows(i, r);
        if (fabs(M(r, lead)) > tol) M.row(r) /= M(r, lead);
        for (int i = 0; i < rowCount; ++i) {
            if (i != r) M.row(i) -= M(i, lead) * M.row(r);
        }
        lead++;
    }
}


//template <unsigned int dim>
//class DirichletAssemblyElasticity : public AssemblyBase<dim>
//{
//public:
//    typedef typename Elasticity::EqFields EqFields;
//    typedef typename Elasticity::EqData EqData;
//
//    static constexpr const char * name() { return "DirichletAssemblyElasticity"; }
//
//    /// Constructor.
//    DirichletAssemblyElasticity(EqFields *eq_fields, EqData *eq_data)
//    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data), drop_tol(1e-12) {
//        this->active_integrals_ = ActiveIntegrals::boundary;
//        this->used_fields_ += eq_fields_->bc_type;
//        this->used_fields_ += eq_fields_->bc_displacement;
//    }
//
//    /// Destructor.
//    ~DirichletAssemblyElasticity() {}
//
//    /// Initialize auxiliary vectors and other data members
//    void initialize(ElementCacheMap *element_cache_map) {
//        this->element_cache_map_ = element_cache_map;
//
//        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
//        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
//        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
//        fe_values_bdr_side_.initialize(*this->quad_low_, *fe_,
//                update_values | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
//        n_dofs_ = fe_->n_dofs();
//        dof_indices_.resize(n_dofs_);
//        vec_view_bdr_ = &fe_values_bdr_side_.vector_view(0);
//    }
//
//
//    /// Assembles boundary integral.
//    void boundary_side_integral(DHCellSide cell_side)
//    {
//    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
//
//        Side side = cell_side.side();
//        const DHCellAccessor &dh_cell = cell_side.cell();
//        dof_indices_ = dh_cell.get_loc_dof_indices();
//        fe_values_bdr_side_.reinit(side);
//
//        unsigned int bc_type;
//        {
//            auto p_side = *( this->boundary_points(cell_side).begin() );
//    	    auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
//    	    bc_type = eq_fields_->bc_type(p_bdr);
//        }
//
//        unsigned int k = 0;
//        if (bc_type == EqFields::bc_type_displacement)
//        {
//            arma::mat A(n_dofs_,n_dofs_+1);
//            for (auto p : this->boundary_points(cell_side) )
//            {
//                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
//                auto bc_displacement = eq_fields_->bc_displacement(p_bdr);
//                for (unsigned int i=0; i<n_dofs_; i++) {
//                    for (unsigned int j=0; j<n_dofs_; j++)
//                        A(i,j) += arma::dot(vec_view_bdr_->value(j,k), vec_view_bdr_->value(i,k)) *
//                            fe_values_bdr_side_.JxW(k);
//                    A(i,n_dofs_) += arma::dot(bc_displacement, vec_view_bdr_->value(i,k)) *
//                            fe_values_bdr_side_.JxW(k);
//                }
//                ++k;
//            }
//            rref(A, drop_tol);
//            for (unsigned int row_idx=0; row_idx<A.n_rows; row_idx++) {
//                arma::uvec indices = arma::find(arma::abs(A(row_idx, arma::span(0,A.n_cols-2))) > drop_tol);
//                ASSERT_LE(indices.size(), 1).error("Dirichlet constraint applies to more than one dof!");
//                if (indices.size() == 1) {
//                    if (std::find(eq_data_->dirichlet_dofs.begin(), eq_data_->dirichlet_dofs.end(), dof_indices_[indices(0)]) == eq_data_->dirichlet_dofs.end()) {
//                        eq_data_->dirichlet_dofs.push_back(dof_indices_[indices(0)]);
//                        eq_data_->dirichlet_coefs.push_back(A(row_idx,indices(0)));
//                        eq_data_->dirichlet_values.push_back(A(row_idx,A.n_cols-1));
//                        eq_data_->dirichlet_row_starts.push_back( eq_data_->dirichlet_row_starts.back() + indices.size() );
//                    }
//                }
//            }
//
//        }
//        else if (bc_type == EqFields::bc_type_displacement_normal)
//        {
//            arma::mat A(n_dofs_,n_dofs_+1);
//            for (auto p : this->boundary_points(cell_side) )
//            {
//                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
//                normal_vector_ = fe_values_bdr_side_.normal_vector(k);
//                auto bc_displacement = eq_fields_->bc_displacement(p_bdr);
//                for (unsigned int i=0; i<n_dofs_; i++) {
//                    for (unsigned int j=0; j<n_dofs_; j++)
//                        A(i,j) += arma::dot(vec_view_bdr_->value(j,k), normal_vector_) *
//                            arma::dot(vec_view_bdr_->value(i,k), normal_vector_) *
//                            fe_values_bdr_side_.JxW(k);
//                    A(i,n_dofs_) += arma::dot(bc_displacement, normal_vector_) *
//                            arma::dot(vec_view_bdr_->value(i,k),normal_vector_) *
//                            fe_values_bdr_side_.JxW(k);
//                }
//                ++k;
//            }
//            rref(A, drop_tol);
//            for (unsigned int row_idx=0; row_idx<A.n_rows; row_idx++) {
//                arma::uvec indices = arma::find(arma::abs(A(row_idx, arma::span(0,A.n_cols-2))) > drop_tol);
//                // ASSERT_LE(indices.size(), 1).error("Dirichlet constraint applies to more than one dof!");
//                if (indices.size() > 0) {
//                    for (unsigned int i=0; i<indices.size(); i++) {
//                        eq_data_->dirichlet_dofs.push_back(dof_indices_[indices(i)]);
//                        eq_data_->dirichlet_coefs.push_back(A(row_idx,indices(i)));
//                    }
//                    eq_data_->dirichlet_values.push_back(A(row_idx,A.n_cols-1));
//                    eq_data_->dirichlet_row_starts.push_back( eq_data_->dirichlet_row_starts.back() + indices.size() );
//                }
//            }
//        }
//
//    }
//
//
//
//
//
//private:
//    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
//
//    /// Data objects shared with Elasticity
//    EqFields *eq_fields_;
//    EqData *eq_data_;
//
//    /// Sub field set contains fields used in calculation.
//    FieldSet used_fields_;
//
//    unsigned int n_dofs_;                                     ///< Number of dofs
//    FEValues<3> fe_values_bdr_side_;                          ///< FEValues of side (boundary integral) object
//
//    LocDofVec dof_indices_;                             ///< Vector of global DOF indices
//    const FEValuesViews::Vector<3> * vec_view_bdr_;           ///< Vector view in boundary calculation.
//
//    arma::vec3 normal_vector_;
//
//    const double drop_tol;                              // Tolerance to drop small values.
//
//    template < template<IntDim...> class DimAssembly>
//    friend class GenericAssembly;
//
//};



#endif /* ASSEMBLY_ELASTICITY_HH_ */
