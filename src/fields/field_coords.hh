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
 * @file    field_coords.hh
 * @brief
 */

#ifndef FIELD_COORDS_HH_
#define FIELD_COORDS_HH_

#include "fields/field_common.hh"                      // for FieldCommon::T...
#include "fields/field_value_cache.hh"                 // for FieldValueCache

namespace IT=Input::Type;


class FieldCoords : public FieldCommon {
public:

    /// Constructor
    FieldCoords()
      : value_cache_( FieldValueCache<double>(3, 1) )
    {
        this->multifield_ = false;
    }

    IT::Instance get_input_type() override {
        ASSERT(false).error("This method can't be used for FieldCoords");

        IT::Abstract abstract = IT::Abstract();
        IT::Instance inst = IT::Instance( abstract, std::vector<IT::TypeBase::ParameterPair>() );
        return inst;
    }

    IT::Array get_multifield_input_type() override {
        ASSERT(false).error("This method can't be used for FieldCoords");

        IT::Array arr = IT::Array( IT::Integer() );
        return arr;
    }

    void set_mesh(const Mesh &mesh) override {
        this->mesh_ = &mesh;
    }

    bool is_constant(FMT_UNUSED Region reg) override {
        return false;
    }

    bool set_time(FMT_UNUSED const TimeStep &time, FMT_UNUSED LimitSide limit_side) override {
        return false;
    }

    void copy_from(FMT_UNUSED const FieldCommon & other) override {
        ASSERT(false).error("Forbidden method for FieldCoords!")
    }

    void field_output(FMT_UNUSED std::shared_ptr<OutputTime> stream) override {
        ASSERT(false).error("Forbidden method for FieldCoords!")
    }

    void observe_output(FMT_UNUSED std::shared_ptr<Observe> observe) override {
        ASSERT(false).error("Forbidden method for FieldCoords!")
    }

    FieldResult field_result(FMT_UNUSED RegionSet region_set) const override {
        return result_none;
    }

    std::string get_value_attribute() const override {
        double limit = std::numeric_limits<double>::max();
        return fmt::format("{{ \"shape\": [ 3, 1 ], \"type\": \"Integer\", \"limit\": [ {}, {} ] }}", -limit, +limit);
    }

    void set_input_list(FMT_UNUSED const Input::Array &list, FMT_UNUSED const TimeGovernor &tg) override
    {}

    /// Implements FieldCommon::cache_allocate
    void cache_reallocate(const ElementCacheMap &cache_map) override {
        unsigned int new_size = ElementCacheMap::n_cached_elements * cache_map.eval_points()->max_size();
        if (new_size > value_cache_.size()) {
            value_cache_.reinit(new_size);
            value_cache_.resize(new_size);
        }
    }

    /// Implements FieldCommon::cache_update
    void cache_update(ElementCacheMap &cache_map, unsigned int i_reg) override {
        // Modify code of FieldSet::update_coords_caches method.
    }

    /// Implements FieldCommon::set_dependency().
    std::vector<const FieldCommon *> set_dependency(FMT_UNUSED FieldSet &field_set, FMT_UNUSED unsigned int i_reg) override {
        return std::vector<const FieldCommon *>();
    }

private:
    /**
     * Field value data cache
     *
     * See implementation of Field<spacedim, Value>::value_cache_
     */
    FieldValueCache<double> value_cache_;

    const Mesh *mesh_;                 ///< Pointer to the mesh.
};

#endif /* FIELD_COORDS_HH_ */
