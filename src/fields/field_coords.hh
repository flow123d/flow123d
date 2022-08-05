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
#include "fields/eval_points.hh"                       // for EvalPoints
#include "fem/mapping_p1.hh"
#include "mesh/ref_element.hh"
#include "mesh/point.hh"                               // for Point

template <int spacedim> class ElementAccessor;

namespace IT=Input::Type;


/**
 * Specialized field represents coordinate variables ('x', 'y', 'z') of FieldFormula.
 */
class FieldCoords : public FieldCommon {
public:
	typedef typename Space<3>::Point Point;

    /// Constructor
    FieldCoords()
      : value_cache_( FieldValueCache<double>(3, 1) )
    {
        this->multifield_ = false;
    	unsigned int cache_size = 1.1 * CacheMapElementNumber::get();
    	value_cache_.reinit(cache_size);
    	value_cache_.resize(cache_size);
    	this->set_shape(3, 1);
    }

    IT::Instance get_input_type() override {
        ASSERT_PERMANENT(false).error("This method can't be used for FieldCoords");

        IT::Abstract abstract = IT::Abstract();
        IT::Instance inst = IT::Instance( abstract, std::vector<IT::TypeBase::ParameterPair>() );
        return inst;
    }

    IT::Array get_multifield_input_type() override {
        ASSERT_PERMANENT(false).error("This method can't be used for FieldCoords");

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
        ASSERT_PERMANENT(false).error("Forbidden method for FieldCoords!");
    }

    void field_output(FMT_UNUSED std::shared_ptr<OutputTime> stream, FMT_UNUSED OutputTime::DiscreteSpace type) override {
        ASSERT_PERMANENT(false).error("Forbidden method for FieldCoords!");
    }

    FieldResult field_result(FMT_UNUSED RegionSet region_set) const override {
        return result_none;
    }

    std::string get_value_attribute() const override {
        double limit = std::numeric_limits<double>::max();
        return fmt::format("{{ \"shape\": [ 3, 1 ], \"type\": \"Double\", \"limit\": [ {}, {} ] }}", -limit, +limit);
    }

    void set_input_list(FMT_UNUSED const Input::Array &list, FMT_UNUSED const TimeGovernor &tg) override
    {}

    /// Implements FieldCommon::cache_allocate
    void cache_reallocate(FMT_UNUSED const ElementCacheMap &cache_map, FMT_UNUSED unsigned int region_idx) const override
    {}

    /// Implements FieldCommon::cache_update
    void cache_update(ElementCacheMap &cache_map, unsigned int region_patch_idx) const override {
    	std::shared_ptr<EvalPoints> eval_points = cache_map.eval_points();
        unsigned int reg_chunk_begin = cache_map.region_chunk_begin(region_patch_idx);
        unsigned int reg_chunk_end = cache_map.region_chunk_end(region_patch_idx);
        unsigned int region_idx = cache_map.eval_point_data(reg_chunk_begin).i_reg_;
        unsigned int last_element_idx = -1;
        ElementAccessor<3> elm;
    	arma::vec3 coords;
        unsigned int dim = 0;

        const MeshBase *mesh; // Holds bulk or boundary mesh by region_idx
        if (region_idx%2 == 1) mesh = shared_->mesh_;
        else mesh = shared_->mesh_->bc_mesh();

        for (unsigned int i_data = reg_chunk_begin; i_data < reg_chunk_end; ++i_data) { // i_eval_point_data
            unsigned int elm_idx = cache_map.eval_point_data(i_data).i_element_;
            if (elm_idx != last_element_idx) {
                elm = mesh->element_accessor( elm_idx );
                dim = elm.dim();
                last_element_idx = elm_idx;
            }

            unsigned int i_point = cache_map.eval_point_data(i_data).i_eval_point_;
            switch (dim) {
            case 0:
                coords = *elm.node(0);
                break;
            case 1:
                coords = MappingP1<1,3>::project_unit_to_real(RefElement<1>::local_to_bary(eval_points->local_point<1>(i_point)),
                        MappingP1<1,3>::element_map(elm));
                break;
            case 2:
                coords = MappingP1<2,3>::project_unit_to_real(RefElement<2>::local_to_bary(eval_points->local_point<2>(i_point)),
                        MappingP1<2,3>::element_map(elm));
                break;
            case 3:
                coords = MappingP1<3,3>::project_unit_to_real(RefElement<3>::local_to_bary(eval_points->local_point<3>(i_point)),
                        MappingP1<3,3>::element_map(elm));
                break;
            default:
            	coords = arma::vec3("0 0 0"); //Should not happen
            }
            value_cache_.set(i_data) = coords;
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
        return std::vector<const FieldCommon *>();
    }

    /// Returns one value of coordinates in one given point @p.
    inline arma::vec3 const & value(const Point &p, FMT_UNUSED const ElementAccessor<3> &elm) const
    {
        return p;
    }

    inline void value_list(const Armor::array &point_list, const ElementAccessor<3> &elm,
                       std::vector<arma::vec3> &value_list) const
    {
        ASSERT(point_list.n_rows() == 3 && point_list.n_cols() == 1).error("Invalid point size.\n");
        ASSERT_EQ(point_list.size(), value_list.size()).error("Different size of point list and value list.\n");

        for (uint i=0; i<point_list.size(); ++i)
            value_list[i] = this->value(point_list.template vec<3>(i), elm);
    }

    /// Return item of @p value_cache_ given by i_cache_point.
    arma::vec3 operator[] (unsigned int i_cache_point) const
    {
        return this->value_cache_.template vec<3>(i_cache_point);
    }

private:
    /**
     * Field value data cache
     *
     * See implementation of Field<spacedim, Value>::value_cache_
     */
    mutable FieldValueCache<double> value_cache_;
};

#endif /* FIELD_COORDS_HH_ */
