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
 * @file    assembly_observe.hh
 * @brief
 */

#ifndef ASSEMBLY_OBSERVE_HH_
#define ASSEMBLY_OBSERVE_HH_

#include <unordered_map>

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "fem/dofhandler.hh"
#include "fem/element_cache_map.hh"
#include "io/observe.hh"


/**
 * @brief Generic class of observe output assemblation.
 *
 * Class
 *  - holds assemblation structures (EvalPoints, Integral objects, Integral data table).
 *  - associates assemblation objects specified by dimension
 *  - provides general assemble method
 *  - provides methods that allow construction of element patches
 */
template < template<IntDim...> class DimAssembly>
class GenericAssemblyObserve : public GenericAssemblyBase
{
public:
    /// Constructor
    GenericAssemblyObserve( typename DimAssembly<1>::EqFields *eq_fields, const std::unordered_set<string> &observe_fields_list,
            std::shared_ptr<Observe> observe)
    : GenericAssemblyBase(),
      multidim_assembly_(eq_fields, observe_fields_list, observe.get(), &this->asm_internals_), observe_(observe), bulk_integral_data_(20, 10)
    {
        multidim_assembly_[1_d]->create_observe_integrals(integrals_);
        multidim_assembly_[2_d]->create_observe_integrals(integrals_);
        multidim_assembly_[3_d]->create_observe_integrals(integrals_);
        this->asm_internals_.element_cache_map_.init(this->asm_internals_.eval_points_);
        multidim_assembly_[1_d]->initialize();
        multidim_assembly_[2_d]->initialize();
        multidim_assembly_[3_d]->initialize();
    }

    /// Getter to set of assembly objects
    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

    /**
     * @brief General assemble methods.
     *
     * Loops through local cells and calls assemble methods of assembly
     * object of each cells over space dimension.
     */
    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) override {
        START_TIMER( DimAssembly<1>::name() );

        unsigned int i_ep, subset_begin, subset_idx;
        auto &patch_point_data = observe_->patch_point_data();
        for(auto & p_data : patch_point_data) {
            subset_idx = integrals_.bulk_[p_data.i_quad]->get_subset_idx();
        	subset_begin = this->asm_internals_.eval_points_->subset_begin(p_data.i_quad+1, subset_idx);
            i_ep = subset_begin + p_data.i_quad_point;
            DHCellAccessor dh_cell = dh->cell_accessor_from_element(p_data.elem_idx);
            bulk_integral_data_.emplace_back(dh_cell, p_data.i_quad_point);
            this->asm_internals_.element_cache_map_.eval_point_data_.emplace_back(p_data.i_reg, p_data.elem_idx, i_ep, 0);
        }
        bulk_integral_data_.make_permanent();
        this->asm_internals_.element_cache_map_.make_paermanent_eval_points();

        this->reallocate_cache();
        this->asm_internals_.element_cache_map_.create_patch();
        multidim_assembly_[1_d]->eq_fields_->cache_update(this->asm_internals_.element_cache_map_);

        multidim_assembly_[1_d]->assemble_cell_integrals(bulk_integral_data_);
        multidim_assembly_[2_d]->assemble_cell_integrals(bulk_integral_data_);
        multidim_assembly_[3_d]->assemble_cell_integrals(bulk_integral_data_);
        bulk_integral_data_.reset();
        this->asm_internals_.element_cache_map_.clear_element_eval_points_map();
        END_TIMER( DimAssembly<1>::name() );
    }


private:
    /// Calls cache_reallocate method on set of used fields
    inline void reallocate_cache() {
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(this->asm_internals_.element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        // DebugOut() << "Order of evaluated fields (" << DimAssembly<1>::name() << "):" << multidim_assembly_[1_d]->eq_fields_->print_dependency();
    }

    MixedPtr<DimAssembly, 1> multidim_assembly_;                  ///< Assembly object
    std::shared_ptr<Observe> observe_;                            ///< Shared Observe object.
    RevertableList<BulkIntegralData> bulk_integral_data_;         ///< Holds data for computing bulk integrals.
};


template <unsigned int dim>
class AssemblyObserveOutput : public AssemblyBase<dim>
{
public:
    typedef EquationOutput EqFields;

    static constexpr const char * name() { return "AssemblyObserveOutput"; }

    /// Constructor.
    AssemblyObserveOutput(EqFields *eq_fields, const std::unordered_set<string> &observe_fields_list, Observe *observe, AssemblyInternals *asm_internals)
    : AssemblyBase<dim>(), eq_fields_(eq_fields), observe_(observe) {
        this->asm_internals_ = asm_internals;
        offsets_.resize(1.1 * CacheMapElementNumber::get());

        for (auto observe_field : observe_fields_list) {
            auto found_field = eq_fields_->field(observe_field);
            used_fields_ += *found_field;
        }
    }

    /// Destructor.
    ~AssemblyObserveOutput() {}

    /// Initialize auxiliary vectors and other data members
    void initialize() {}

    /// Assembles the cell integrals for the given dimension.
    inline void assemble_cell_integrals(const RevertableList<GenericAssemblyBase::BulkIntegralData> &bulk_integral_data) {
        unsigned int element_patch_idx, field_value_cache_position, val_idx;
        this->reset_offsets();
        for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
            if (bulk_integral_data[i].cell.dim() != dim) continue;
            element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(bulk_integral_data[i].cell.elm_idx());
            auto p = *( bulk_integral_->points(element_patch_idx).begin()); // evaluation point
            field_value_cache_position = this->asm_internals_->element_cache_map_.element_eval_point(element_patch_idx, p.eval_point_idx() + bulk_integral_data[i].subset_index);
            val_idx = ObservePointAccessor(observe_, i).loc_point_time_index();
            this->offsets_[field_value_cache_position] = val_idx;
        }
        for (FieldListAccessor f_acc : this->used_fields_.fields_range()) {
            f_acc->fill_observe_value(observe_->get_output_cache(f_acc->name()), this->offsets_);
        }
    }


    /// Create bulk integral according to dim
    void create_observe_integrals(AssemblyIntegrals &integrals) {
        std::vector<arma::vec> reg_points;

        auto &patch_point_data = observe_->patch_point_data();
        for(auto & p_data : patch_point_data) {
            auto el_acc = eq_fields_->mesh()->element_accessor(p_data.elem_idx);
            if (el_acc.dim()!=dim) continue;
            p_data.i_reg = el_acc.region_idx().idx();
            p_data.i_quad = el_acc.dim() - 1;
            p_data.i_quad_point = reg_points.size();
            reg_points.push_back(p_data.local_coords);
        }

        if (reg_points.size() > 0) {
            this->quad_ = new Quadrature(dim, reg_points.size());
            for (uint j=0; j<reg_points.size(); j++) {
                arma::vec::fixed<dim> fix_p = reg_points[j].subvec(0, dim-1);
                this->quad_->weight(j) = 1.0;
                this->quad_->set(j) = fix_p;
            }
            bulk_integral_ = this->create_bulk_integral(this->quad_);
            this->integrals_.bulk_ = bulk_integral_;
            integrals.bulk_[dim-1] = bulk_integral_;
        }
    }

private:
    void reset_offsets() {
        std::fill(offsets_.begin(), offsets_.end(), -1);
    }

    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    Observe *observe_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output
    std::vector<int> offsets_;                                ///< Holds indices (offsets) of cached data to output data vector

    std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;     ///< Accessor of integral

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssemblyObserve;
};




#endif /* ASSEMBLY_OBSERVE_HH_ */
