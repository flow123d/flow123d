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
 * @file    field_depth.hh
 * @brief
 */

#ifndef FIELD_DEPTH_HH_
#define FIELD_DEPTH_HH_

#include "fields/field_common.hh"                      // for FieldCommon::T...
#include "fields/field_value_cache.hh"                 // for FieldValueCache
#include "fields/eval_points.hh"                       // for EvalPoints
#include "fields/field_coords.hh"
#include "fields/surface_depth.hh"
#include "fem/mapping_p1.hh"
#include "mesh/ref_element.hh"

namespace IT=Input::Type;


/**
 * Specialized field represents surface depth ('d' variable) of FieldFormula.
 */
class FieldDepth : public FieldCommon {
public:

    /// Constructor
	FieldDepth()
      : value_cache_( FieldValueCache<double>(1, 1) ), surface_depth_(nullptr)
    {
        this->multifield_ = false;
    	unsigned int cache_size = 1.1 * CacheMapElementNumber::get();
    	value_cache_.reinit(cache_size);
    	value_cache_.resize(cache_size);
    	this->set_shape(1, 1);
    }

    IT::Instance get_input_type() override {
        ASSERT(false).error("This method can't be used for FieldDepth");

        IT::Abstract abstract = IT::Abstract();
        IT::Instance inst = IT::Instance( abstract, std::vector<IT::TypeBase::ParameterPair>() );
        return inst;
    }

    IT::Array get_multifield_input_type() override {
        ASSERT(false).error("This method can't be used for FieldDepth");

        IT::Array arr = IT::Array( IT::Integer() );
        return arr;
    }

    void set_mesh(const Mesh &mesh) override {
        shared_->mesh_ = &mesh;
    }

    bool is_constant(FMT_UNUSED Region reg) override {
        return false;
    }

    bool set_time(FMT_UNUSED const TimeStep &time, FMT_UNUSED LimitSide limit_side) override {
        return false;
    }

    void copy_from(FMT_UNUSED const FieldCommon & other) override {
        ASSERT(false).error("Forbidden method for FieldCoords!");
    }

    void field_output(FMT_UNUSED std::shared_ptr<OutputTime> stream, FMT_UNUSED OutputTime::DiscreteSpace type) override {
        ASSERT(false).error("Forbidden method for FieldCoords!");
    }

    void observe_output(FMT_UNUSED std::shared_ptr<Observe> observe) override {
        ASSERT(false).error("Forbidden method for FieldCoords!");
    }

    FieldResult field_result(FMT_UNUSED RegionSet region_set) const override {
        return result_none;
    }

    std::string get_value_attribute() const override {
        double limit = std::numeric_limits<double>::max();
        return fmt::format("{{ \"shape\": [ 1, 1 ], \"type\": \"Double\", \"limit\": [ {}, {} ] }}", -limit, +limit);
    }

    void set_input_list(FMT_UNUSED const Input::Array &list, FMT_UNUSED const TimeGovernor &tg) override
    {}

    /// Implements FieldCommon::cache_allocate
    void cache_reallocate(FMT_UNUSED const ElementCacheMap &cache_map, FMT_UNUSED unsigned int region_idx) const override
    {}

    /// Implements FieldCommon::cache_update
    void cache_update(ElementCacheMap &cache_map, unsigned int region_patch_idx) const override {
    	if (surface_depth_ == nullptr) return;

    	std::shared_ptr<EvalPoints> eval_points = cache_map.eval_points();
        unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_patch_idx);
        unsigned int reg_chunk_end = cache_map.region_chunk_end(region_patch_idx);
        auto * coords_cache = field_coords_->value_cache();
    	arma::vec3 p; // evaluated point

        for (unsigned int i_data = reg_chunk_begin; i_data < reg_chunk_end; ++i_data) { // i_eval_point_data
            p = coords_cache->template vec<3>(i_data);
            value_cache_.set(i_data) = surface_depth_->compute_distance(p);
        }
    }

    /// Implements FieldCommon::value_cache
    FieldValueCache<double> * value_cache() override {
    	return &value_cache_;
    }

    /// Implements FieldCommon::value_cache
    const FieldValueCache<double> * value_cache() const override {
    	return &value_cache_;
    }

    /// Implements FieldCommon::set_dependency().
    std::vector<const FieldCommon *> set_dependency(FMT_UNUSED FieldSet &field_set, FMT_UNUSED unsigned int i_reg) const override {
    	std::vector<const FieldCommon *> res;
    	res.push_back(field_coords_);
        return res;
    }

    /// Setter of surface_depth data member
    inline void set_surface_depth(std::shared_ptr<SurfaceDepth> surface_depth) {
        surface_depth_ = surface_depth;
    }

    /// Setter of field_coords data member
    inline void set_field_coords(FieldCoords * field_coords) {
    	field_coords_ = field_coords;
    }

    /// Forbidden method in this class.
    OutputTime::OutputDataPtr output_data_cache(FMT_UNUSED OutputTime::DiscreteSpace space_type, FMT_UNUSED std::shared_ptr<OutputTime> stream) const override {
    	ASSERT(false);
    	return nullptr;
    }

    /// Forbidden method in this class.
    void fill_data_value(FMT_UNUSED BulkPoint &p, FMT_UNUSED unsigned int value_idx,
                         FMT_UNUSED std::shared_ptr<ElementDataCacheBase> output_data_base) override
    {
    	ASSERT(false);
    }

private:
    /**
     * Field value data cache
     *
     * See implementation of Field<spacedim, Value>::value_cache_
     */
    mutable FieldValueCache<double> value_cache_;

    /// Surface depth object calculate distance from surface.
    std::shared_ptr<SurfaceDepth> surface_depth_;

    FieldCoords * field_coords_;        ///< Pointer to coordinates field.
};

#endif /* FIELD_DEPTH_HH_ */
