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

