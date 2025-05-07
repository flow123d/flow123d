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
#include "fem/element_cache_map.hh"
#include "io/output_time.hh"
#include "io/element_data_cache.hh"
#include "mesh/ref_element.hh"



/**
 * Common ancestor of AssemblyOutputElemData and AssemblyOutputNodeData class.
 */
template <unsigned int dim>
class AssemblyOutputBase : public AssemblyBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;

    /// Constructor.
    AssemblyOutputBase(unsigned int quad_order, EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(quad_order), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        offsets_.resize(CacheMapElementNumber::get());
    }

    /// Constructor.
    AssemblyOutputBase(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        offsets_.resize(CacheMapElementNumber::get());
    }

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }

    /// Sets output data members.
    void set_output_data(const FieldSet &used, std::shared_ptr<OutputTime> stream) {
    	used_fields_ = FieldSet();
    	used_fields_ += used;
    	stream_ = stream;
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        for (FieldListAccessor f_acc : used_fields_.fields_range()) {
            stream_->update_time(f_acc->time());
        }
    }


protected:
    void reset_offsets() {
        std::fill(offsets_.begin(), offsets_.end(), -1);
    }

    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    EqData *eq_data_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output
    std::shared_ptr<OutputTime> stream_;                      ///< Output stream.
    std::vector<int> offsets_;                                ///< Holds indices (offsets) of cached data to output data vector
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyOutputElemData : public AssemblyOutputBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;
    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;

    static constexpr const char * name() { return "AssemblyOutputElemData"; }

    /// Constructor.
    AssemblyOutputElemData(EqFields *eq_fields, EqData *eq_data)
    : AssemblyOutputBase<dim>(0, eq_fields, eq_data) {}

    /// Destructor.
    ~AssemblyOutputElemData() {}

    /// Assembles the cell integrals for the given dimension.
    inline void assemble_cell_integrals() {
    	if (dim!=1) return;  // Perform full output in one loop
    	unsigned int element_patch_idx, field_value_cache_position, val_idx;
    	this->reset_offsets();
    	for (unsigned int i=0; i<this->bulk_integral_data_.permanent_size(); ++i) {
            element_patch_idx = this->element_cache_map_->position_in_cache(this->bulk_integral_data_[i].cell.elm_idx());
            auto p = *( this->bulk_points(element_patch_idx).begin() ); // evaluation point (in element center)
            field_value_cache_position = this->element_cache_map_->element_eval_point(element_patch_idx, p.eval_point_idx());
            val_idx = this->stream_->get_output_mesh_ptr()->get_loc_elem_idx(this->bulk_integral_data_[i].cell.elm_idx());
            this->offsets_[field_value_cache_position] = val_idx;
    	}
        for (FieldListAccessor f_acc : this->used_fields_.fields_range()) {
            f_acc->fill_data_value(this->offsets_);
        }
    }

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyOutputNodeData : public AssemblyOutputBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;
    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;

    static constexpr const char * name() { return "AssemblyOutputNodeData"; }

    /// Constructor.
    AssemblyOutputNodeData(EqFields *eq_fields, EqData *eq_data)
    : AssemblyOutputBase<dim>(eq_fields, eq_data) {
        this->quad_ = new Quadrature(dim, RefElement<dim>::n_nodes);
        for(unsigned int i = 0; i<RefElement<dim>::n_nodes; i++)
        {
            this->quad_->weight(i) = 1.0;
            this->quad_->set(i) = RefElement<dim>::node_coords(i);
        }
    }

    /// Destructor.
    ~AssemblyOutputNodeData() {}

    /// Sets output data members.
    void set_output_data(const FieldSet &used, std::shared_ptr<OutputTime> stream) {
    	AssemblyOutputBase<dim>::set_output_data(used, stream);
    	offset_vec_ = this->stream_->get_output_mesh_ptr()->offsets()->get_data();
    }


    /// Assembles the cell integrals for the given dimension.
    inline void assemble_cell_integrals() {
    	if (dim!=1) return;  // Perform full output in one loop
    	unsigned int element_patch_idx, field_value_cache_position, val_idx;
    	this->reset_offsets();
    	for (unsigned int i=0; i<this->bulk_integral_data_.permanent_size(); ++i) {
            element_patch_idx = this->element_cache_map_->position_in_cache(this->bulk_integral_data_[i].cell.elm_idx());
            val_idx = (*offset_vec_)[ this->stream_->get_output_mesh_ptr()->get_loc_elem_idx(this->bulk_integral_data_[i].cell.elm_idx()) ];
            auto p = *( this->bulk_points(element_patch_idx).begin() );
            field_value_cache_position = this->element_cache_map_->element_eval_point(element_patch_idx, p.eval_point_idx());
            for (uint j=0; j<this->bulk_integral_data_[i].cell.dim()+1; ++j) {
                this->offsets_[field_value_cache_position+j] = val_idx+j;
            }
    	}
        for (FieldListAccessor f_acc : this->used_fields_.fields_range()) {
            f_acc->fill_data_value(this->offsets_);
        }
    }

private:
    std::shared_ptr< std::vector<unsigned int> > offset_vec_;   ///< Holds offsets of individual local elements

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

#endif /* ASSEMBLY_OUTPUT_HH_ */
