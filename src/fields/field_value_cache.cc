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
#include "fields/eval_subset.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"


/******************************************************************************
 * Implementation of FieldValueCache methods
 */

template<class elm_type, class return_type>
FieldValueCache<elm_type, return_type>::FieldValueCache(unsigned int n_rows, unsigned int n_cols)
: data_(0, n_rows, n_cols), eval_points_(nullptr), dim_(EvalPoints::undefined_dim) {
    used_subsets_.fill(false);
    subset_starts_.fill(-1);
}

template<class elm_type, class return_type>
FieldValueCache<elm_type, return_type>::~FieldValueCache() {}

template<class elm_type, class return_type>
void FieldValueCache<elm_type, return_type>::init(std::shared_ptr<EvalSubset> eval_subset, unsigned int n_cache_elements) {
    ASSERT_EQ(dim_, EvalPoints::undefined_dim).error("Repeated initialization!\n");

    this->n_cache_elements_ = n_cache_elements;
    eval_points_ = eval_subset->eval_points();
    data_.resize(n_cache_elements * eval_points_->size());
    subset_starts_[0] = 0;
    for (uint i=0; i<eval_points_->n_subsets(); ++i)
    	subset_starts_[i+1] = eval_points_->subset_end(i) * n_cache_elements;
    dim_ = eval_points_->point_dim();
}

template<class elm_type, class return_type>
void FieldValueCache<elm_type, return_type>::mark_used(std::shared_ptr<EvalSubset> sub_set) {
    used_subsets_[sub_set->get_subset_idx()] = true;
}


/******************************************************************************
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::n_cached_elements = 20;

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();


ElementCacheMap::ElementCacheMap()
: elm_idx_(ElementCacheMap::n_cached_elements, ElementCacheMap::undef_elem_idx),
  dim_(EvalPoints::undefined_dim), ready_to_reading_(false) {}


void ElementCacheMap::add(const DHCellAccessor &dh_cell) {
	ASSERT_DBG(ready_to_reading_);
    ASSERT_LT(update_data_.added_elements_.size(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    unsigned int elm_idx = dh_cell.elm_idx();
    update_data_.added_elements_.insert(elm_idx);
}


void ElementCacheMap::prepare_elements_to_update(Mesh *mesh) {
	// Start of cache update
	ready_to_reading_ = false;

	// Find elements that were stored in previous cache update.
    std::array<bool, ElementCacheMap::n_cached_elements> stored_previous;
	unsigned int n_stored_element = 0;
    stored_previous.fill(false);
    for (auto it = update_data_.added_elements_.begin(); it != update_data_.added_elements_.end();)
        if (cache_idx_.find(*it) != cache_idx_.end()) {
           stored_previous[cache_idx_.at(*it)] = true;
           ++n_stored_element;
           it = update_data_.added_elements_.erase(it); // erase element from added_elements set that exists in cache
        } else it++;

    // Prepare preserved data for move to beginning part of data block, update in elm_idx_, cache_idx_
    for (unsigned int i_source=n_stored_element, i_target=0; i_source<stored_previous.size(); ++i_source) {
        if (elm_idx_[i_source] == ElementCacheMap::undef_elem_idx) break;
        if (!stored_previous[i_source]) {
            cache_idx_.erase( cache_idx_.find(elm_idx_[i_source]) );
            elm_idx_[i_source] = ElementCacheMap::undef_elem_idx;
        	continue;
        }
        while (stored_previous[i_target] && (i_target<n_stored_element)) i_target++;
        update_data_.preserved_elements_[i_source] = i_target;
        cache_idx_[ elm_idx_[i_source] ] = i_target;
        cache_idx_.erase( cache_idx_.find(elm_idx_[i_source]) );
        elm_idx_[i_target] = elm_idx_[i_source];
        elm_idx_[i_source] = ElementCacheMap::undef_elem_idx;
        i_target++;
    }

    // Distribute elements by region
    for (const auto &elm_idx : update_data_.added_elements_) {
        ElementAccessor<3> elm(mesh, elm_idx);
        unsigned int reg_idx = elm.region_idx().idx();
        typename std::unordered_map<unsigned int, ElementSet>::iterator region_it = update_data_.region_element_map_.find(reg_idx);
        if (region_it == update_data_.region_element_map_.end()) {
        	update_data_.region_element_map_.insert( {reg_idx, ElementSet()} );
            region_it = update_data_.region_element_map_.find(reg_idx);
        }
        region_it->second.push_back( elm );
    }

    // Set new elements to elm_idx_, cache_idx_ sorted by region
    ASSERT_EQ_DBG(cache_idx_.size(), n_stored_element);
    for (auto region_it = update_data_.region_element_map_.begin(); region_it != update_data_.region_element_map_.end(); region_it++) {
        update_data_.region_cache_begin_[region_it->first] = cache_idx_.size();
        for (auto elm : region_it->second) {
            unsigned int elm_idx = elm.idx();
            cache_idx_[elm_idx] = n_stored_element;
            elm_idx_[n_stored_element] = elm_idx;
            n_stored_element++;
        }
    }
}


void ElementCacheMap::clear_elements_to_update() {
	update_data_.added_elements_.clear();
	update_data_.preserved_elements_.clear();
	update_data_.region_element_map_.clear();
	update_data_.region_cache_begin_.clear();
	ready_to_reading_ = true; // end of cache update
}


DHCellAccessor & ElementCacheMap::operator() (DHCellAccessor &dh_cell) const {
	ASSERT_DBG(ready_to_reading_);
	unsigned int elm_idx = dh_cell.elm_idx();
	std::unordered_map<unsigned int, unsigned int>::const_iterator it = cache_idx_.find(elm_idx);
	if ( it != cache_idx_.end() ) dh_cell.set_element_cache_index( it->second );
	else dh_cell.set_element_cache_index( ElementCacheMap::undef_elem_idx );
    return dh_cell;
}


/******************************************************************************
 * Explicit instantiation of templates
 */

template class FieldValueCache<unsigned int, unsigned int>;    // NOT tested, necessary for linking!
template class FieldValueCache<int, int>;             // NOT tested, necessary for linking!
template class FieldValueCache<double, double>;
template class FieldValueCache<double, arma::vec3>;
template class FieldValueCache<double, arma::mat33>;
