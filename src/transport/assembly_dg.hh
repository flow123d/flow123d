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
 * @file    assembly_dg.hh
 * @brief
 */

#ifndef ASSEMBLY_DG_HH_
#define ASSEMBLY_DG_HH_

#include "transport/advection_diffusion_model.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"


/**
 * Base (non-templated) class of container for Finite element and related objects.
 */
class AssemblyDGBase
{
public:
    typedef std::vector<std::shared_ptr<AssemblyDGBase> > MultidimAssemblyDG;

    virtual ~AssemblyDGBase() {}

    virtual void initialize() = 0;

    virtual void assemble_mass_matrix(DHCellAccessor cell, std::vector<Vec> &ret_vec, LinSys **ls_dt) = 0;
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyDG : public AssemblyDGBase
{
public:

    /// Constructor.
    AssemblyDG(unsigned int fe_order, AdvectionDiffusionModel &adm)
    : fe_(new FE_P_disc<dim>(fe_order)), fe_rt_(new FE_RT0<dim>),
	  quad_(new QGauss<dim>(2*fe_order)), mapping_(new MappingP1<dim,3>),
	  model_(adm), fe_values_(*mapping_, *quad_, *fe_, update_values | update_JxW_values | update_quadrature_points) {

    	ndofs_ = fe_->n_dofs();
        qsize_ = quad_->size();
        dof_indices_.resize(ndofs_);
    }

    /// Destructor.
    ~AssemblyDG() {
        delete fe_;
        delete fe_rt_;
        delete quad_;
        delete mapping_;
    }

    /// Initialize auxiliary vectors and other data members
    void initialize() override {
        local_mass_matrix_.resize(ndofs_*ndofs_);
        local_retardation_balance_vector_.resize(ndofs_);
        local_mass_balance_vector_.resize(ndofs_);

        mm_coef_.resize(qsize_);
        ret_coef_.resize(model_.n_substances());
        for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
        {
            ret_coef_[sbi].resize(qsize_);
        }
    }

    /// Assemble integral over element
    void assemble_mass_matrix(DHCellAccessor cell, std::vector<Vec> &ret_vec, LinSys **ls_dt) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");
        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        model_.compute_mass_matrix_coefficient(fe_values_.point_list(), elm, mm_coef_);
        model_.compute_retardation_coefficient(fe_values_.point_list(), elm, ret_coef_);

        for (unsigned int sbi=0; sbi<model_.n_substances(); ++sbi)
        {
            // assemble the local mass matrix
            for (unsigned int i=0; i<ndofs_; i++)
            {
                for (unsigned int j=0; j<ndofs_; j++)
                {
                    local_mass_matrix_[i*ndofs_+j] = 0;
                    for (unsigned int k=0; k<qsize_; k++)
                        local_mass_matrix_[i*ndofs_+j] += (mm_coef_[k]+ret_coef_[sbi][k])*fe_values_.shape_value(j,k)*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                }
            }

            for (unsigned int i=0; i<ndofs_; i++)
            {
                local_mass_balance_vector_[i] = 0;
                local_retardation_balance_vector_[i] = 0;
                for (unsigned int k=0; k<qsize_; k++)
                {
                    local_mass_balance_vector_[i] += mm_coef_[k]*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                    local_retardation_balance_vector_[i] -= ret_coef_[sbi][k]*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                }
            }

            model_.balance()->add_mass_matrix_values(model_.get_subst_idx()[sbi], elm.region().bulk_idx(), dof_indices_, local_mass_balance_vector_);
            ls_dt[sbi]->mat_set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_mass_matrix_[0]));
            VecSetValues(ret_vec[sbi], ndofs_, &(dof_indices_[0]), &(local_retardation_balance_vector_[0]), ADD_VALUES);
        }
    }

private:

    FiniteElement<dim> *fe_;     ///< Finite element for the solution of the advection-diffusion equation.
    FiniteElement<dim> *fe_rt_;  ///< Finite element for the water velocity field.
    Quadrature<dim> *quad_;      ///< Quadrature used in assembling methods.
    MappingP1<dim,3> *mapping_;  ///< Auxiliary mapping of reference elements.

    /// Reference to model (we must use common ancestor of concentration and heat model)
    AdvectionDiffusionModel &model_;

    unsigned int ndofs_;                                      ///< Number of dofs
    unsigned int qsize_;                                      ///< Size of FEValues quadrature
    FEValues<dim,3> fe_values_;                               ///< FEValues of object (assemble_mass_matrix method)
    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<PetscScalar> local_mass_matrix_;                   ///< Helper vector for assembly mass matrix
    vector<PetscScalar> local_retardation_balance_vector_;    ///< Same as previous.
    vector<PetscScalar> local_mass_balance_vector_;           ///< Same as previous.

	/// Mass matrix coefficients.
	vector<double> mm_coef_;
	/// Retardation coefficient due to sorption.
	vector<vector<double> > ret_coef_;

    friend class FEObjects;
};



#endif /* FE_VALUE_HANDLER_HH_ */
