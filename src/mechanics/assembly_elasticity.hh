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
    typedef typename Elasticity::EqData EqDataDG;

    /// Constructor.
    StiffnessAssemblyElasticity(EqDataDG *data)
    : AssemblyBase<dim>(1), data_(data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
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
        qsize_ = this->quad_->size();
        qsize_low_ = this->quad_low_->size();
        dof_indices_.resize(n_dofs_);
        local_matrix_.resize(n_dofs_*n_dofs_);
        vec_view_ = &fe_values_.vector_view(0);
        vec_view_side_ = &fe_values_side_.vector_view(0);
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
        if (bc_type == EqDataDG::bc_type_displacement)
        {
            unsigned int k=0;
            for (auto p : data_->stiffness_assembly_->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += data_->dirichlet_penalty(p) * side_measure *
                                arma::dot(vec_view_side_->value(i,k),vec_view_side_->value(j,k)) * fe_values_side_.JxW(k);
                k++;
            }
        }
        else if (bc_type == EqDataDG::bc_type_displacement_normal)
        {
            unsigned int k=0;
            for (auto p : data_->stiffness_assembly_->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += data_->dirichlet_penalty(p) * side_measure *
                                arma::dot(vec_view_side_->value(i,k), fe_values_side_.normal_vector(k)) *
                                arma::dot(vec_view_side_->value(j,k), fe_values_side_.normal_vector(k)) * fe_values_side_.JxW(k);
                k++;
            }
        }

        data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }


    /// Assembles the fluxes between elements of different dimensions.
    inline void neigbour_integral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        // MOVE Elasticity::assemble_matrix_element_side code

    	/* Used data members:
    	 * fe_values_side_, fe_values_sub_, n_dofs_ (replace ndofs_side), n_dofs_sub_, qsize_low_
        vector<vector<int> > side_dof_indices(2);
        vector<unsigned int> n_dofs = { ndofs_sub, ndofs_side };
    	vector<double> frac_sigma(qsize);
    	vector<double> csection_lower(qsize), csection_higher(qsize), young(qsize), poisson(qsize), alpha(qsize);
        PetscScalar local_matrix[2][2][(ndofs_side)*(ndofs_side)];
    	*/

    	/*if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        // Note: use data members csection_ and velocity_ for appropriate quantities of lower dim element

        double comm_flux[2][2];
        unsigned int n_dofs[2];
        ElementAccessor<3> elm_lower_dim = cell_lower_dim.elm();
        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i] = dof_indices_[i];
        }
        fe_values_vb_.reinit(elm_lower_dim);
        n_dofs[0] = fv_sb_[0]->n_dofs();

        DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i+n_dofs[0]] = dof_indices_[i];
        }
        fe_values_side_.reinit(neighb_side.side());
        n_dofs[1] = fv_sb_[1]->n_dofs();

        */
    }


    /// Implements @p AssemblyBase::reallocate_cache.
    void reallocate_cache(const ElementCacheMap &cache_map) override
    {
        data_->cache_reallocate(cache_map);
    }


    private:
        shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
        shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

        /// Data object shared with EqDataDG
        EqDataDG *data_;

        unsigned int n_dofs_;                                     ///< Number of dofs
        unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
        unsigned int qsize_;                                      ///< Size of quadrature
        unsigned int qsize_low_;                                  ///< Size of quadrature of dim-1
        FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
        FEValues<3> fe_values_side_;                              ///< FEValues of side object
        FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
        const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
        const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in boundary side integral calculation.

        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};



#endif /* ASSEMBLY_ELASTICITY_HH_ */

