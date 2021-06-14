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
 * @file    assembly_convection.hh
 * @brief
 */

#ifndef ASSEMBLY_CONVECTION_HH_
#define ASSEMBLY_CONVECTION_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "transport/transport.h"
//#include "fem/fe_p.hh"
//#include "fem/fe_values.hh"
//#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class MassAssemblyConvection : public AssemblyBase<dim>
{
public:
    typedef typename ConvectionTransport::EqFields EqFields;
    typedef typename ConvectionTransport::EqData EqData;

    static constexpr const char * name() { return "MassAssemblyConvection"; }

    /// Constructor.
    MassAssemblyConvection(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->water_content;
    }

    /// Destructor.
    ~MassAssemblyConvection() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elm = cell.elm();
        // we have currently zero order P_Disc FE
        ASSERT_DBG(cell.get_loc_dof_indices().size() == 1);
        IntIdx local_p0_dof = cell.get_loc_dof_indices()[0];

        auto p = *( this->bulk_points(element_patch_idx).begin() );
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
            eq_data_->balance_->add_mass_values(eq_data_->subst_idx[sbi], cell, {local_p0_dof},
                    {eq_fields_->cross_section(p)*eq_fields_->water_content(p)*elm.measure()}, 0);

        VecSetValue(eq_data_->mass_diag, eq_data_->dh_->get_local_to_global_map()[local_p0_dof],
                eq_fields_->cross_section(p)*eq_fields_->water_content(p), INSERT_VALUES);
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        VecZeroEntries(eq_data_->mass_diag);
        eq_data_->balance_->start_mass_assembly(eq_data_->subst_idx);
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        eq_data_->balance_->finish_mass_assembly(eq_data_->subst_idx);

        VecAssemblyBegin(eq_data_->mass_diag);
        VecAssemblyEnd(eq_data_->mass_diag);

        eq_data_->is_mass_diag_changed = true;
    }

    private:
        shared_ptr<FiniteElement<dim>> fe_;                    ///< Finite element for the solution of the advection-diffusion equation.

        /// Data objects shared with TransportDG
        EqFields *eq_fields_;
        EqData *eq_data_;

        /// Sub field set contains fields used in calculation.
        FieldSet used_fields_;

        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};


#endif /* ASSEMBLY_CONVECTION_HH_ */
