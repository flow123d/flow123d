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
#include "fields/field_value_cache.hh"
#include "io/observe.hh"


/// Holds data of one eval point on patch (index of element and local coordinations).
struct PatchPointData {
    /// Default constructor
    PatchPointData() {}

    /// Constructor with data mebers initialization
    PatchPointData(unsigned int elm_idx, arma::vec loc_coords)
    : elem_idx(elm_idx), local_coords(loc_coords), i_quad(0), i_quad_point(0) {}

    /// Copy constructor
    PatchPointData(const PatchPointData &other)
    : elem_idx(other.elem_idx), local_coords(other.local_coords),
      i_quad(other.i_quad), i_quad_point(other.i_quad_point) {}

    unsigned int elem_idx;        ///< Index of element
    arma::vec local_coords;       ///< Local coords of point
    unsigned int i_reg;           ///< Index of region (use during patch creating)
    unsigned int i_quad;          ///< Index of quadrature (use during patch creating), i_quad = dim-1
    unsigned int i_quad_point;    ///< Index of point in quadrature (use during patch creating)
};
typedef std::vector<PatchPointData> PatchPointVec;


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
    : multidim_assembly_(eq_fields, eq_data), bulk_integral_data_(20, 10)
    {
        eval_points_ = std::make_shared<EvalPoints>();
        multidim_assembly_[1_d]->create_observe_integrals(eval_points_, integrals_);
        multidim_assembly_[2_d]->create_observe_integrals(eval_points_, integrals_);
        multidim_assembly_[3_d]->create_observe_integrals(eval_points_, integrals_);
        element_cache_map_.init(eval_points_);
        multidim_assembly_[1_d]->initialize(&element_cache_map_);
        multidim_assembly_[2_d]->initialize(&element_cache_map_);
        multidim_assembly_[3_d]->initialize(&element_cache_map_);
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
        auto &patch_point_data = multidim_assembly_[1_d]->eq_data_->patch_point_data();
        for(auto & p_data : patch_point_data) {
            subset_idx = integrals_.bulk_[p_data.i_quad]->get_subset_idx();
        	subset_begin = eval_points_->subset_begin(p_data.i_quad+1, subset_idx);
            i_ep = subset_begin + p_data.i_quad_point;
            DHCellAccessor dh_cell = dh->cell_accessor_from_element(p_data.elem_idx);
            bulk_integral_data_.emplace_back(dh_cell, p_data.i_quad_point);
            element_cache_map_.eval_point_data_.emplace_back(p_data.i_reg, p_data.elem_idx, i_ep, 0);
        }
        bulk_integral_data_.make_permanent();
        element_cache_map_.eval_point_data_.make_permanent();

        this->reallocate_cache();
        element_cache_map_.create_patch();
        multidim_assembly_[1_d]->eq_fields_->cache_update(element_cache_map_);

        multidim_assembly_[1_d]->assemble_cell_integrals(bulk_integral_data_);
        multidim_assembly_[2_d]->assemble_cell_integrals(bulk_integral_data_);
        multidim_assembly_[3_d]->assemble_cell_integrals(bulk_integral_data_);
        bulk_integral_data_.reset();
        END_TIMER( DimAssembly<1>::name() );
    }


private:
    /// Calls cache_reallocate method on set of used fields
    inline void reallocate_cache() {
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(this->element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        // DebugOut() << "Order of evaluated fields (" << DimAssembly<1>::name() << "):" << multidim_assembly_[1_d]->eq_fields_->print_dependency();
    }

    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    MixedPtr<DimAssembly, 1> multidim_assembly_;                  ///< Assembly object
    AssemblyIntegrals integrals_;                                 ///< Holds integral objects.
    RevertableList<BulkIntegralData> bulk_integral_data_;         ///< Holds data for computing bulk integrals.
};


template <unsigned int dim>
class AssemblyObserveOutput : public AssemblyBase<dim>
{
public:
    typedef EquationOutput EqFields;
    typedef EquationOutput EqData;

    static constexpr const char * name() { return "AssemblyObserveOutput"; }

    /// Constructor.
    AssemblyObserveOutput(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(), eq_fields_(eq_fields), eq_data_(eq_data) {
        offsets_.resize(1.1 * CacheMapElementNumber::get());

        for (auto observe_field : eq_data_->observe_fields()) {
            auto found_field = eq_fields_->field(observe_field);
            used_fields_ += *found_field;
        }

        if (dim==1) {
            auto observe_ptr = eq_data_->output_stream()->observe( eq_data_->eq_mesh() );
            auto &patch_point_data = eq_data_->patch_point_data();
            for (ObservePointAccessor op_acc : observe_ptr->local_range()) {
                patch_point_data.emplace_back(op_acc.observe_point().element_idx(), op_acc.observe_point().local_coords());
            }
        }
    }

    /// Destructor.
    ~AssemblyObserveOutput() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }

    /// Assembles the cell integrals for the given dimension.
    inline void assemble_cell_integrals(const RevertableList<GenericAssemblyBase::BulkIntegralData> &bulk_integral_data) {
        unsigned int element_patch_idx, field_value_cache_position, val_idx;
        this->reset_offsets();
        for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
            if (bulk_integral_data[i].cell.dim() != dim) continue;
            element_patch_idx = this->element_cache_map_->position_in_cache(bulk_integral_data[i].cell.elm_idx());
            auto p = *( this->bulk_points(element_patch_idx).begin()); // evaluation point
            field_value_cache_position = this->element_cache_map_->element_eval_point(element_patch_idx, p.eval_point_idx() + bulk_integral_data[i].subset_index);
            val_idx = ObservePointAccessor(eq_data_->output_stream()->observe( eq_data_->eq_mesh() ).get(), i).loc_point_time_index();
            this->offsets_[field_value_cache_position] = val_idx;
        }
        for (FieldListAccessor f_acc : this->used_fields_.fields_range()) {
            f_acc->fill_observe_value(this->offsets_);
        }
    }


    /// Create bulk integral according to dim
    void create_observe_integrals(std::shared_ptr<EvalPoints> eval_points, AssemblyIntegrals &integrals) {
        std::vector<arma::vec> reg_points;

        auto &patch_point_data = eq_data_->patch_point_data();
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
            this->integrals_.bulk_ = eval_points->add_bulk<dim>(*this->quad_);
            integrals.bulk_[dim-1] = this->integrals_.bulk_;
        }
    }

private:
    void reset_offsets() {
        std::fill(offsets_.begin(), offsets_.end(), -1);
    }

    /// Data objects shared with EquationOutput
    EqFields *eq_fields_;
    EqData *eq_data_;

    FieldSet used_fields_;                                    ///< Sub field set contains fields performed to output
    std::vector<int> offsets_;                                ///< Holds indices (offsets) of cached data to output data vector


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssemblyObserve;
};




#endif /* ASSEMBLY_OBSERVE_HH_ */
