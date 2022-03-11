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
//#include "fem/fe_values.hh"
//#include "quadrature/quadrature_lib.hh"
//#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"
#include "io/observe.hh"
//#include "io/element_data_cache.hh"
//#include "mesh/ref_element.hh"

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
	GenericAssemblyObserve( typename DimAssembly<1>::EqFields *eq_fields, typename DimAssembly<1>::EqData *eq_data)
    : multidim_assembly_(eq_fields, eq_data), observe_(nullptr)
    {
        eval_points_ = std::make_shared<EvalPoints>();
        element_cache_map_.init(eval_points_);
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
        ASSERT_PTR(observe_).error("Uninitialized observe object!\n");

        START_TIMER( DimAssembly<1>::name() );
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(this->element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        this->create_patch(dh);
        multidim_assembly_[1_d]->eq_fields_->cache_update(element_cache_map_);

        multidim_assembly_[1_d]->assemble_cell_integrals();
        END_TIMER( DimAssembly<1>::name() );
    }

    inline void set_observe(std::shared_ptr<Observe> observe) {
        this->observe_ = observe;
    }


private:

    void create_patch(std::shared_ptr<DOFHandlerMultiDim> dh) {
        eval_points_->clear();
        element_cache_map_.eval_point_data_.reset();

        unsigned int reg_idx, ele_idx, i_ep;
        MeshBase *dh_mesh = dh->mesh();
        for(ObservePointAccessor op_acc : observe_->local_range()) {
        	ele_idx = op_acc.observe_point().element_idx();
            auto el_acc = dh_mesh->element_accessor(ele_idx);
            reg_idx = el_acc.region_idx().idx();
            switch (el_acc.dim()) {
            case 1:
                i_ep = eval_points_->add_local_point<1>(op_acc.observe_point().local_coords());
                break;
            case 2:
                i_ep = eval_points_->add_local_point<2>(op_acc.observe_point().local_coords());
                break;
            case 3:
                i_ep = eval_points_->add_local_point<3>(op_acc.observe_point().local_coords());
                break;
            default:
                ASSERT_PERMANENT(false);
                break;
            }
            element_cache_map_.eval_point_data_.emplace_back(reg_idx, ele_idx, i_ep, op_acc.loc_point_time_index());
            // Value dh_loc_idx_ holds index of observe point in data cache.
        }

        element_cache_map_.eval_point_data_.make_permanent();
        element_cache_map_.create_patch();
    }

    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    MixedPtr<DimAssembly, 1> multidim_assembly_;                  ///< Assembly object
    std::shared_ptr<Observe> observe_;                            ///< Shared observe object
};


template <unsigned int dim>
class AssemblyObserveOutput : public AssemblyBase<dim>
{
public:
    // TODO set used_fields: same as in AssemblyOutputBase
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;
    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;

    static constexpr const char * name() { return "AssemblyObserveOutput"; }

    /// Constructor.
    AssemblyObserveOutput(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(), eq_fields_(eq_fields), eq_data_(eq_data) {}

    /// Destructor.
    ~AssemblyObserveOutput() {}

    /// Sets observe output data members (set of used fields).
    void set_observe_data(const FieldSet &used) {
    	used_fields_ = FieldSet();
    	used_fields_ += used;
    }

    /// Assembles the cell integrals for the given dimension.
    inline void assemble_cell_integrals(/*const RevertableList<BulkIntegralData> &bulk_integral_data*/) {
        if (dim!=1) return;  // Perform full output in one loop
//        unsigned int element_patch_idx, field_value_cache_position, val_idx;
//        this->reset_offsets();
//        for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
//            element_patch_idx = this->element_cache_map_->position_in_cache(bulk_integral_data[i].cell.elm_idx());
//            auto p = *( this->bulk_points(element_patch_idx).begin() ); // evaluation point (in element center)
//            field_value_cache_position = this->element_cache_map_->element_eval_point(element_patch_idx, p.eval_point_idx());
//            val_idx = this->stream_->get_output_mesh_ptr()->get_loc_elem_idx(bulk_integral_data[i].cell.elm_idx());
//            this->offsets_[field_value_cache_position] = val_idx;
//        }
//        for (FieldListAccessor f_acc : this->used_fields_.fields_range()) {
//            f_acc->fill_data_value(this->offsets_);
//        }
    }

private:
    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    EqData *eq_data_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssemblyObserve;
};




#endif /* ASSEMBLY_OBSERVE_HH_ */
