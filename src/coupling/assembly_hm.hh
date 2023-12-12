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
 * @file    assembly_hm.hh
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
    typedef typename HM_Iterative::EqData EqData;
    typedef typename DarcyLMH::EqFields FlowEqFields;

    static constexpr const char * name() { return "FlowPotentialAssemblyHM"; }

    /// Constructor.
    FlowPotentialAssemblyHM(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->alpha;
        this->used_fields_ += eq_fields_->density;
        this->used_fields_ += eq_fields_->gravity;
        this->used_fields_ += eq_data_->flow_->eq_fields().bc_type;
        this->used_fields_ += eq_data_->flow_->eq_fields().bc_pressure;
    }

    /// Destructor.
    ~FlowPotentialAssemblyHM() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        shared_ptr<FiniteElement<dim>> fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_values_side_.initialize(*this->quad_low_, *fe_, update_side_JxW_values);

        dof_indices_.resize(fe_->n_dofs());

        ref_potential_vec_ = eq_fields_->ref_potential_ptr_->vec();
    }


    /// Assemble integral over element
    inline void boundary_side_integral(DHCellSide dh_side)
    {
        dof_indices_ = dh_side.cell().get_loc_dof_indices();

        double ref_pot = 0;
        if (dh_side.side().is_boundary())
        {
            auto p_side = *this->boundary_points(dh_side).begin();
            auto p_bdr = p_side.point_bdr( dh_side.cond().element_accessor() );
            unsigned int flow_bc_type = eq_data_->flow_->eq_fields().bc_type(p_bdr);
            if (flow_bc_type == DarcyLMH::EqFields::dirichlet || flow_bc_type == DarcyLMH::EqFields::total_flux)
            {
                unsigned int k=0;
                for ( auto p : this->boundary_points(dh_side) )
                {
                    // The reference potential is applied only on dirichlet and total_flux b.c.,
                    // i.e. where only mechanical traction is prescribed.
                    auto p_bdr = p.point_bdr(dh_side.cond().element_accessor());
                    double alpha = eq_fields_->alpha(p);
                    double density = eq_fields_->density(p);
                    double gravity = eq_fields_->gravity(p);
                    double bc_pressure = eq_data_->flow_->eq_fields().bc_pressure(p_bdr);
                    ref_pot += -alpha*density*gravity*bc_pressure * fe_values_side_.JxW(k) / dh_side.measure();
                    ++k;
                }
            }
        }
        ref_potential_vec_.set(dof_indices_[dh_side.side_idx()], ref_pot);
    }



private:

    EqFields *eq_fields_;       ///< Fields shared with HM_Iterative
    EqData *eq_data_;           ///< Data objects shared with HM_Iterative

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    FEValues<3> fe_values_side_;                              ///< FEValues of side object
    LocDofVec dof_indices_;                             ///< Vector of global DOF indices
    VectorMPI ref_potential_vec_;                             ///< Vector of dofs of field ref_potential

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class ResidualAssemblyHM : public AssemblyBase<dim>
{
public:
    typedef typename HM_Iterative::EqFields EqFields;
    typedef typename HM_Iterative::EqData EqData;

    static constexpr const char * name() { return "ResidualAssemblyHM"; }

    /// Constructor.
    ResidualAssemblyHM(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk);
        this->used_fields_ += eq_data_->flow_->eq_fields().field_ele_pressure;
        this->used_fields_ += eq_fields_->old_iter_pressure;
    }

    /// Destructor.
    ~ResidualAssemblyHM() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
        shared_ptr<FE_P<dim>> fe_ = std::make_shared< FE_P<dim> >(0);
        fe_values_.initialize(*this->quad_, *fe_, update_JxW_values);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;

        fe_values_.reinit(cell.elm());

        // compute pressure error
        unsigned int k=0;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            double new_p = eq_data_->flow_->eq_fields().field_ele_pressure(p);
            double old_p = eq_fields_->old_iter_pressure(p);
            eq_data_->p_dif2 += (new_p - old_p)*(new_p - old_p) * fe_values_.JxW(k);
            eq_data_->p_norm2 += old_p*old_p * fe_values_.JxW(k);
            ++k;
        }
    }



private:
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    FEValues<3> fe_values_;                                   ///< FEValues of cell object

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ASSEMBLY_HM_HH_ */

