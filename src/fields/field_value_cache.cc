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
 * @file    field_value_cache.cc
 * @brief
 * @author  David Flanderka
 */

#include <limits>
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/eval_points.hh"
#include "fem/dh_cell_accessor.hh"


/******************************************************************************
 * Implementation of FieldValueCache methods
 */

template<class Value>
FieldValueCache<Value>::FieldValueCache()
: data_(0, Value::NRows_, Value::NCols_) {}

template<class Value>
FieldValueCache<Value>::FieldValueCache(EvalPoints eval_points)
: data_(ElementCacheMap::n_cached_elements * eval_points.size(), Value::NRows_, Value::NCols_) {
	dim_ = eval_points.point_dim();
}

template<class Value>
FieldValueCache<Value>::~FieldValueCache() {}

template<class Value>
void FieldValueCache<Value>::mark_used(EvalSubset sub_quad) {
    auto local_point_vec = sub_quad.get_point_indices(0);
    for (auto p_idx : local_point_vec) used_points_.insert(p_idx);
}


/******************************************************************************
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::n_cached_elements = 20;

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();


ElementCacheMap::ElementCacheMap(unsigned int dim)
: elm_idx_(ElementCacheMap::n_cached_elements, ElementCacheMap::undef_elem_idx),
  begin_idx_(0), end_idx_(0), dim_(dim) {}


unsigned int ElementCacheMap::add(DHCellAccessor dh_cell) {
	ASSERT_LT(added_elements_.size(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is to long!\n");
	unsigned int elm_idx = dh_cell.elm_idx();
	std::unordered_map<unsigned int, unsigned int>::iterator it = cache_idx_.find(elm_idx);
    if ( it != cache_idx_.end() ) {
        return it->second;
    } else {
    	ASSERT( std::find(added_elements_.begin(), added_elements_.end(), elm_idx) == added_elements_.end() )(elm_idx).error("Repeated addition of element!\n");

    	added_elements_.push_back(elm_idx);
    	unsigned int idx_pos = end_idx_;
    	end_idx_ = (end_idx_+1) % ElementCacheMap::n_cached_elements;
        return idx_pos;
    }
}


void ElementCacheMap::clear_elements_to_update() {
	begin_idx_ = end_idx_;
	added_elements_.clear();
}


DHCellAccessor & ElementCacheMap::operator() (DHCellAccessor &dh_cell) const {
	unsigned int elm_idx = dh_cell.elm_idx();
	std::unordered_map<unsigned int, unsigned int>::const_iterator it = cache_idx_.find(elm_idx);
	if ( it != cache_idx_.end() ) dh_cell.set_element_cache_index( it->second );
	else dh_cell.set_element_cache_index( ElementCacheMap::undef_elem_idx );
    return dh_cell;
}


/******************************************************************************
 * Explicit instantiation of templates
 */

template class FieldValueCache<FieldValue<0>::Enum >;    // NOT tested, necessary for linking!
template class FieldValueCache<FieldValue<0>::Integer >; // NOT tested, necessary for linking!
template class FieldValueCache<FieldValue<0>::Scalar >;
template class FieldValueCache<FieldValue<3>::VectorFixed >;
template class FieldValueCache<FieldValue<3>::TensorFixed >;
