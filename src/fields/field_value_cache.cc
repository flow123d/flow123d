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
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();
const unsigned int ElementCacheMap::simd_size_double = 4;


ElementCacheMap::ElementCacheMap()
: elm_idx_(ElementCacheMap::n_cached_elements, ElementCacheMap::undef_elem_idx),
  ready_to_reading_(false), element_eval_points_map_(nullptr), eval_point_data_(0),
  regions_starts_(2*ElementCacheMap::regions_in_chunk,ElementCacheMap::regions_in_chunk),
  element_starts_(2*ElementCacheMap::elements_in_chunk,ElementCacheMap::elements_in_chunk) {
    cache_idx_.reserve(ElementCacheMap::n_cached_elements);
}


ElementCacheMap::~ElementCacheMap() {
    if (element_eval_points_map_!=nullptr) {
        delete element_eval_points_map_;
    }
}


void ElementCacheMap::init(std::shared_ptr<EvalPoints> eval_points) {
    this->eval_points_ = eval_points;
    unsigned int ep_data_size = ElementCacheMap::n_cached_elements * eval_points_->max_size();
    eval_point_data_.resize(ep_data_size);
    element_eval_points_map_ = new int [ElementCacheMap::n_cached_elements * eval_points->max_size()];
}


void ElementCacheMap::add(const DHCellAccessor &dh_cell) {
    /// obsolete method
	ASSERT_DBG(!ready_to_reading_);
    ASSERT_LT_DBG(this->n_elements(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    this->add_to_region(dh_cell.elm());
}


void ElementCacheMap::add(const DHCellSide &cell_side) {
    /// obsolete method
	ASSERT_DBG(!ready_to_reading_);
    ASSERT_LT_DBG(this->n_elements(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    this->add_to_region(cell_side.cell().elm());
}


void ElementCacheMap::add(const ElementAccessor<3> &elm_acc) {
    /// obsolete method
	ASSERT_DBG(!ready_to_reading_);
    ASSERT_LT_DBG(this->n_elements(), ElementCacheMap::n_cached_elements).error("ElementCacheMap overflowed. List of added elements is too long!\n");
    this->add_to_region(elm_acc);
}


void ElementCacheMap::prepare_elements_to_update() {
    std::sort(eval_point_data_.begin(), eval_point_data_.end());
    unsigned int last_region_idx = -1;
    unsigned int last_element_idx = -1;
    unsigned int i_pos=0; // position in eval_point_data_
    regions_starts_.reset();
    element_starts_.reset();
    regions_to_map_.clear();
    element_to_map_.clear();
    for (auto it=eval_point_data_.begin(); it!=eval_point_data_.end(); ++it, ++i_pos) {
        if (it->i_element_ != last_element_idx) { // new element
            if (it->i_reg_ != last_region_idx) { // new region
                regions_to_map_[it->i_reg_] = regions_starts_.temporary_size();
                regions_starts_.push_back( element_starts_.temporary_size() );
                last_region_idx = it->i_reg_;
            }
            element_to_map_[it->i_element_] = element_starts_.temporary_size();
            element_starts_.push_back(i_pos);
            last_element_idx = it->i_element_;
        }
    }
    regions_starts_.push_back( element_starts_.temporary_size() );
    element_starts_.push_back(i_pos);
    regions_starts_.make_permanent();
    element_starts_.make_permanent();

    /*** OLD CODE of create map ***/
    // Erase element data of previous step
    cache_idx_.clear();
    std::fill(elm_idx_.begin(), elm_idx_.end(), ElementCacheMap::undef_elem_idx);
    this->clear_element_eval_points_map();

    // Set new elements to elm_idx_, cache_idx_ sorted by region
	unsigned int n_stored_element = 0, n_region = 0;
	update_data_.region_element_cache_range_[0] = 0;
    for (auto region_it = update_data_.region_cache_indices_map_.begin(); region_it != update_data_.region_cache_indices_map_.end(); region_it++) {
    	region_it->second.cache_position_ = n_region;
        for (unsigned int i_elm=0; i_elm<region_it->second.n_elements_; ++i_elm) {
            unsigned int elm_idx = region_it->second.elm_indices_[i_elm];
            cache_idx_[elm_idx] = n_stored_element;
            elm_idx_[n_stored_element] = elm_idx;
            n_stored_element++;
        }
        update_data_.region_element_cache_range_[n_region+1] = n_stored_element;
        n_region++;
    }
    /*** end of OLD CODE ***/
}


void ElementCacheMap::create_elements_points_map() {
    /*** OLD CODE of create map ***/
    unsigned int size = this->eval_points_->max_size();
    unsigned int idx_to_region = 1;
    unsigned int region_last_elm = update_data_.region_element_cache_range_[idx_to_region];
    unsigned int points_in_cache = 0;
    update_data_.region_value_cache_range_[0] = 0;
	for (unsigned int i_elm=0; i_elm<ElementCacheMap::n_cached_elements; ++i_elm) {
	    for (unsigned int i_point=0; i_point<size; ++i_point) {
	        if (element_eval_point(i_elm, i_point) == ElementCacheMap::point_in_proggress) {
	            set_element_eval_point(i_elm, i_point, points_in_cache);
	            points_in_cache++;
	        }
	    }
        if (region_last_elm==i_elm+1) {
            while (points_in_cache%ElementCacheMap::simd_size_double > 0) points_in_cache++;
            update_data_.region_value_cache_range_[idx_to_region] = points_in_cache;
            idx_to_region++;
            region_last_elm = update_data_.region_element_cache_range_[idx_to_region];
        }
	}
    /*** end of OLD CODE ***/
}


void ElementCacheMap::start_elements_update() {
	ready_to_reading_ = false;
}

void ElementCacheMap::finish_elements_update() {
	update_data_.region_cache_indices_map_.clear(); /*** OLD CODE of create map ***/
	ready_to_reading_ = true;
}

void ElementCacheMap::mark_used_eval_points(const DHCellAccessor &dh_cell, unsigned int subset_idx, unsigned int data_size, unsigned int start_point) {
    /// obsolete method
    unsigned int elem_idx_in_cache = cache_idx_[dh_cell.elm().mesh_idx()];
    unsigned int points_begin = eval_points_->subset_begin(dh_cell.dim(), subset_idx) + start_point;
    for (unsigned int i=points_begin; i<points_begin+data_size; ++i)
        set_element_eval_point(elem_idx_in_cache, i, ElementCacheMap::point_in_proggress);
}


void ElementCacheMap::mark_used_eval_points(const ElementAccessor<3> elm, unsigned int subset_idx, unsigned int data_size, unsigned int start_point) {
    /// obsolete method
    unsigned int elem_idx_in_cache = cache_idx_[elm.mesh_idx()];
    unsigned int points_begin = eval_points_->subset_begin(elm.dim(), subset_idx) + start_point;
    for (unsigned int i=points_begin; i<points_begin+data_size; ++i)
        set_element_eval_point(elem_idx_in_cache, i, ElementCacheMap::point_in_proggress);
}


void ElementCacheMap::clear_element_eval_points_map() {
    /// obsolete method
	ASSERT_PTR_DBG(element_eval_points_map_);
    unsigned int size = this->eval_points_->max_size();
	for (unsigned int i_elm=0; i_elm<ElementCacheMap::n_cached_elements; ++i_elm)
	    for (unsigned int i_point=0; i_point<size; ++i_point)
	        set_element_eval_point(i_elm, i_point, ElementCacheMap::unused_point);
}


void ElementCacheMap::add_to_region(ElementAccessor<3> elm) {
    /// obsolete method

    unsigned int reg_idx = elm.region_idx().idx();
    typename std::unordered_map<unsigned int, RegionData>::iterator region_it = update_data_.region_cache_indices_map_.find(reg_idx);
    if (region_it == update_data_.region_cache_indices_map_.end()) {
    	update_data_.region_cache_indices_map_.insert( {reg_idx, RegionData()} );
        region_it = update_data_.region_cache_indices_map_.find(reg_idx);
    }

    region_it->second.add(elm);
}


DHCellAccessor & ElementCacheMap::cache_map_index(DHCellAccessor &dh_cell) const {
	ASSERT_DBG(ready_to_reading_);
	unsigned int elm_idx = dh_cell.elm_idx();
	std::unordered_map<unsigned int, unsigned int>::const_iterator it = cache_idx_.find(elm_idx);
	if ( it != cache_idx_.end() ) dh_cell.set_element_cache_index( it->second );
	else dh_cell.set_element_cache_index( ElementCacheMap::undef_elem_idx );
    return dh_cell;
}

