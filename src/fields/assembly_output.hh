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

#ifndef ASSEMBLY_OUTPUT_HH_
#define ASSEMBLY_OUTPUT_HH_

#include <unordered_map>

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"
#include "io/output_time.hh"
#include "io/element_data_cache.hh"
#include "mesh/ref_element.hh"



typedef typename std::unordered_map<std::string, OutputTime::OutputDataPtr> UsedElementDataCaches;


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyOutputElemData : public AssemblyBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;

    static constexpr const char * name() { return "AssemblyOutputElemData"; }

    /// Constructor.
    AssemblyOutputElemData(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
    }

    /// Destructor.
    ~AssemblyOutputElemData() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }

    /// Sets output data members.
    void set_output_data(const FieldSet &used, UsedElementDataCaches used_caches, std::shared_ptr<OutputTime> stream) {
    	used_fields_ = FieldSet();
    	used_fields_ += used;
    	used_caches_ = used_caches;
    	stream_ = stream;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        auto p = *( this->bulk_points(element_patch_idx).begin() ); // evaluation point (in element center)
        unsigned int val_idx = stream_->get_output_mesh_ptr()->get_loc_elem_idx(cell.elm_idx());
        for (FieldListAccessor f_acc : used_fields_.fields_range()) {
            typename OutputTime::OutputDataPtr output_data_base = used_caches_[f_acc->name()];
            f_acc->fill_data_value(p, val_idx, output_data_base);
        }
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        for (FieldListAccessor f_acc : used_fields_.fields_range()) {
            stream_->update_time(f_acc->time());
        }
    }

private:
    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    EqData *eq_data_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output
    UsedElementDataCaches used_caches_;                       ///< Map of ElementDataCaches assigned to items of used_fields_
    std::shared_ptr<OutputTime> stream_;                      ///< Output stream.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyOutputNodeData : public AssemblyBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;

    static constexpr const char * name() { return "AssemblyOutputNodeData"; }

    /// Constructor.
    AssemblyOutputNodeData(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;

        this->quad_ = new Quadrature(dim, RefElement<dim>::n_nodes);
        for(unsigned int i = 0; i<RefElement<dim>::n_nodes; i++)
        {
            this->quad_->weight(i) = 1.0;
            this->quad_->set(i) = RefElement<dim>::node_coords(i);
        }
    }

    /// Destructor.
    ~AssemblyOutputNodeData() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }

    /// Sets output data members.
    void set_output_data(const FieldSet &used, UsedElementDataCaches used_caches, std::shared_ptr<OutputTime> stream) {
    	used_fields_ = FieldSet();
    	used_fields_ += used;
    	used_caches_ = used_caches;
    	stream_ = stream;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        unsigned int output_elm_idx = stream_->get_output_mesh_ptr()->get_loc_elem_idx(cell.elm_idx());
        auto &offset_vec = *( stream_->get_output_mesh_ptr()->offsets()->get_component_data(0).get() );
        unsigned int offset = offset_vec[output_elm_idx];
        for (auto p : this->bulk_points(element_patch_idx) ) {
            for (FieldListAccessor f_acc : used_fields_.fields_range()) {
                typename OutputTime::OutputDataPtr output_data_base = used_caches_[f_acc->name()];
                f_acc->fill_data_value(p, offset, output_data_base);
            }
            ++offset;
        }
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        for (FieldListAccessor f_acc : used_fields_.fields_range()) {
            stream_->update_time(f_acc->time());
        }
    }

private:
    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    EqData *eq_data_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output
    UsedElementDataCaches used_caches_;                       ///< Map of ElementDataCaches assigned to items of used_fields_
    std::shared_ptr<OutputTime> stream_;                      ///< Output stream.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

#endif /* ASSEMBLY_OUTPUT_HH_ */
