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
: data_(0, Value::NRows_, Value::NCols_), eval_points_(nullptr) {}

template<class Value>
FieldValueCache<Value>::FieldValueCache(const EvalPoints *eval_points, unsigned int n_cache_points)
: data_(n_cache_points * eval_points->size(), Value::NRows_, Value::NCols_),
  eval_points_(eval_points), dim_(eval_points->point_dim()) {}

template<class Value>
FieldValueCache<Value>::~FieldValueCache() {}

template<class Value>
void FieldValueCache<Value>::mark_used(EvalSubset sub_set) {
    for (unsigned int i=0; i<sub_set.n_permutations(); ++i)
	    used_points_.insert(sub_set.get_block_idx(i));
}


/******************************************************************************
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::n_cached_elements = 20;

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();


ElementCacheMap::ElementCacheMap(unsigned int dim)
: elm_idx_(ElementCacheMap::n_cached_elements, ElementCacheMap::undef_elem_idx),
  begin_idx_(0), end_idx_(0), dim_(dim) {}


void ElementCacheMap::add(DHCellAccessor dh_cell) {
    ASSERT_LT(added_elements_.size(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is to long!\n");
    unsigned int elm_idx = dh_cell.elm_idx();
   	added_elements_.insert(elm_idx);
}


void ElementCacheMap::prepare_elements_to_update() {
    // Compute number of element stored in data cache (stored in previous cache update).
    unsigned int n_stored_element = 0;
    for (auto elm_idx : added_elements_)
        if (cache_idx_.find(elm_idx) != cache_idx_.end()) ++n_stored_element;

    // Test if new elements can be add the end of cache
    if (cache_idx_.size() + added_elements_.size() - n_stored_element <= ElementCacheMap::n_cached_elements) {
    	// Erase elements from added_elements set that exist in cache
        for (auto it = added_elements_.begin(); it != added_elements_.end();) {
            if ( cache_idx_.find(*it) != cache_idx_.end() ) it = added_elements_.erase(it);
            else it++;
        }
        begin_idx_ = cache_idx_.size();
    } else {
    	// Clear cache, there is not sufficient space for new elements
        std::fill(elm_idx_.begin(), elm_idx_.end(), ElementCacheMap::undef_elem_idx);
        cache_idx_.clear();
        begin_idx_ = 0;
    }
    end_idx_ = begin_idx_ + added_elements_.size();

    // Add new elements indices to cache_idx_ and elm_idx_
    unsigned int cache_pos = begin_idx_;
    for (auto el_idx : added_elements_) {
    	cache_idx_[el_idx] = cache_pos;
		elm_idx_[cache_pos] = el_idx;
        cache_pos++;
    }
}


void ElementCacheMap::clear_elements_to_update() {
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
