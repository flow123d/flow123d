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
 * @file    assembly_lmh.hh
 * @brief
 */

#ifndef ASSEMBLY_LMH_HH_
#define ASSEMBLY_LMH_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "flow/darcy_flow_lmh.hh"
//#include "fem/fe_p.hh"
//#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


template <unsigned int dim>
class ReadInitCondAssemblyLMH : public AssemblyBase<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyLMH::EqData EqData;

    static constexpr const char * name() { return "ReadInitCondAssemblyLMH"; }

    /// Constructor.
    ReadInitCondAssemblyLMH(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->init_pressure;
    }

    /// Destructor.
    virtual ~ReadInitCondAssemblyLMH() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        p_indices_ = cell.cell_with_other_dh(eq_data_->dh_p_.get()).get_loc_dof_indices();
        ASSERT_DBG(p_indices_.n_elem == 1);
        l_indices_ = cell.cell_with_other_dh(eq_data_->dh_cr_.get()).get_loc_dof_indices();

		// set initial condition
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        init_value_ = eq_fields_->init_pressure(p);
        p_idx_ = eq_data_->dh_p_->parent_indices()[p_indices_[0]];
        eq_data_->full_solution.set(p_idx_, init_value_);

        for (unsigned int i=0; i<cell.elm()->n_sides(); i++) {
             init_value_on_edge_ = init_value_ / cell.elm().side(i)->edge().n_sides();
             l_idx_ = eq_data_->dh_cr_->parent_indices()[l_indices_[i]];
             eq_data_->full_solution.add(l_idx_, init_value_on_edge_);

             eq_data_->p_edge_solution.add(l_indices_[i], init_value_on_edge_);
        }

        update_water_content(cell, p);
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        eq_data_->full_solution.ghost_to_local_begin();
        eq_data_->full_solution.ghost_to_local_end();

        eq_data_->p_edge_solution.ghost_to_local_begin();
        eq_data_->p_edge_solution.ghost_to_local_end();
        eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
    }



protected:
    virtual void update_water_content(FMT_UNUSED const DHCellAccessor& dh_cell, FMT_UNUSED BulkPoint &p) {}

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

private:
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    LocDofVec p_indices_, l_indices_;            ///< Vectors of local DOF indices pre-computed on different DOF handlers
    unsigned int p_idx_, l_idx_;                 ///< Local DOF indices extract from previous vectors
    double init_value_, init_value_on_edge_;     ///< Pre-computed values of init_pressure.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

#endif /* ASSEMBLY_LMH_HH_ */

