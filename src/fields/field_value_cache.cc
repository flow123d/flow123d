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
#include "fields/field_value_cache.impl.hh"
#include "fields/field_values.hh"
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/accessors.hh"


/******************************************************************************
 * Implementation of FieldValueCache methods
 */

template<class elm_type>
FieldValueCache<elm_type>::FieldValueCache(unsigned int n_rows, unsigned int n_cols)
: data_(n_rows, n_cols) {}

template<class elm_type>
FieldValueCache<elm_type>::~FieldValueCache() {}

template<class elm_type>
void FieldValueCache<elm_type>::init(std::shared_ptr<EvalPoints> eval_points, unsigned int n_cache_elements) {
    ASSERT_EQ(data_.size(), 0).error("Repeated initialization!\n");

    this->n_cache_points_ = n_cache_elements * eval_points->max_size();
    data_.reinit(n_cache_points_);
    data_.resize(n_cache_points_);
}


/******************************************************************************
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();


ElementCacheMap::ElementCacheMap()
: elm_idx_(ElementCacheMap::n_cached_elements, ElementCacheMap::undef_elem_idx),
  ready_to_reading_(false), element_eval_points_map_(nullptr), points_in_cache_(0) {
    cache_idx_.reserve(ElementCacheMap::n_cached_elements);
    update_data_.n_elements_ = 0;
}


ElementCacheMap::~ElementCacheMap() {
    if (element_eval_points_map_!=nullptr) {
        for (unsigned int i=0; i<ElementCacheMap::n_cached_elements; ++i)
            delete element_eval_points_map_[i];
        delete element_eval_points_map_;
    }
}


void ElementCacheMap::init(std::shared_ptr<EvalPoints> eval_points) {
	this->eval_points_ = eval_points;

	unsigned int size = this->eval_points_->max_size();
	element_eval_points_map_ = new int* [ElementCacheMap::n_cached_elements];
	for (unsigned int i=0; i<ElementCacheMap::n_cached_elements; ++i)
	    element_eval_points_map_[i] = new int [size];
}


void ElementCacheMap::add(const DHCellAccessor &dh_cell) {
	ASSERT_DBG(!ready_to_reading_);
    ASSERT_LT(update_data_.n_elements_, ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    this->add_to_region(dh_cell.elm());
}


void ElementCacheMap::add(const DHCellSide &cell_side) {
	ASSERT_DBG(!ready_to_reading_);
    ASSERT_LT(update_data_.n_elements_, ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    this->add_to_region(cell_side.cell().elm());
}


void ElementCacheMap::prepare_elements_to_update() {
    // Erase element data of previous step
    cache_idx_.clear();
    std::fill(elm_idx_.begin(), elm_idx_.end(), ElementCacheMap::undef_elem_idx);
    this->clear_element_eval_points_map();

    // Set new elements to elm_idx_, cache_idx_ sorted by region
	unsigned int n_stored_element = 0, n_region = 0;
	update_data_.region_element_cache_range_[0] = 0;
    for (auto region_it = update_data_.region_cache_indices_map_.begin(); region_it != update_data_.region_cache_indices_map_.end(); region_it++) {
    	update_data_.region_cache_indices_range_[region_it->first] = n_region;
        for (unsigned int i_elm=0; i_elm<region_it->second.n_elements_; ++i_elm) {
            unsigned int elm_idx = region_it->second.elm_indices_[i_elm];
            cache_idx_[elm_idx] = n_stored_element;
            elm_idx_[n_stored_element] = elm_idx;
            n_stored_element++;
        }
        update_data_.region_element_cache_range_[n_region+1] = n_stored_element;
        n_region++;
    }
}


void ElementCacheMap::create_elements_points_map() {
    unsigned int size = this->eval_points_->max_size();
    unsigned int idx_to_region = 1;
    unsigned int region_last_elm = update_data_.region_element_cache_range_[idx_to_region];
    points_in_cache_ = 0;
    update_data_.region_value_cache_range_[0] = 0;
	for (unsigned int i_elm=0; i_elm<ElementCacheMap::n_cached_elements; ++i_elm) {
	    for (unsigned int i_point=0; i_point<size; ++i_point) {
	        if (element_eval_points_map_[i_elm][i_point] == ElementCacheMap::point_in_proggress) {
	    	    element_eval_points_map_[i_elm][i_point] = points_in_cache_;
	            points_in_cache_++;
	        }
	    }
        if (region_last_elm==i_elm+1) {
            update_data_.region_value_cache_range_[idx_to_region] = points_in_cache_;
        	idx_to_region++;
        	region_last_elm = update_data_.region_element_cache_range_[idx_to_region];
        }
	}
}


void ElementCacheMap::start_elements_update() {
	update_data_.n_elements_ = 0;
	ready_to_reading_ = false;
}

void ElementCacheMap::finish_elements_update() {
	update_data_.region_cache_indices_map_.clear();
	update_data_.region_cache_indices_range_.clear();
	ready_to_reading_ = true;
}

void ElementCacheMap::mark_used_eval_points(const DHCellAccessor &dh_cell, unsigned int subset_idx, unsigned int data_size, unsigned int start_point) {
    unsigned int elem_idx_in_cache = cache_idx_[dh_cell.elm_idx()];
    unsigned int points_begin = eval_points_->subset_begin(dh_cell.dim(), subset_idx) + start_point;
    for (unsigned int i=points_begin; i<points_begin+data_size; ++i)
        element_eval_points_map_[elem_idx_in_cache][i] = ElementCacheMap::point_in_proggress;
}


void ElementCacheMap::clear_element_eval_points_map() {
	ASSERT_PTR_DBG(element_eval_points_map_);
    unsigned int size = this->eval_points_->max_size();
	for (unsigned int i_elm=0; i_elm<ElementCacheMap::n_cached_elements; ++i_elm)
	    for (unsigned int i_point=0; i_point<size; ++i_point)
	        element_eval_points_map_[i_elm][i_point] = ElementCacheMap::unused_point;
}


void ElementCacheMap::add_to_region(ElementAccessor<3> elm) {
    // TODO: see unordered_map::insert(map::begin(), pair) variant
    // this can be used and avoid dupicit find and the condition.

    unsigned int reg_idx = elm.region_idx().idx();
    typename std::unordered_map<unsigned int, RegionData>::iterator region_it = update_data_.region_cache_indices_map_.find(reg_idx);
    if (region_it == update_data_.region_cache_indices_map_.end()) {
    	update_data_.region_cache_indices_map_.insert( {reg_idx, RegionData()} );
        region_it = update_data_.region_cache_indices_map_.find(reg_idx);
    }

    // TODO: What is the reason for this condition, the elm should not be duplicate.
    if ( region_it->second.add(elm) ) update_data_.n_elements_++;
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

template class FieldValueCache<unsigned int>;
template class FieldValueCache<int>;
template class FieldValueCache<double>;
