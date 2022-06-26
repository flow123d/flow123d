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
#include "mesh/accessors.hh"


/******************************************************************************
 * Implementation of ElementCacheMap methods
 */

const unsigned int ElementCacheMap::undef_elem_idx = std::numeric_limits<unsigned int>::max();
const unsigned int ElementCacheMap::simd_size_double = 4;


ElementCacheMap::ElementCacheMap()
: elm_idx_(CacheMapElementNumber::get(), ElementCacheMap::undef_elem_idx),
  ready_to_reading_(false), element_eval_points_map_(nullptr), eval_point_data_(0),
  regions_starts_(2*ElementCacheMap::regions_in_chunk,ElementCacheMap::regions_in_chunk),
  element_starts_(2*ElementCacheMap::elements_in_chunk,ElementCacheMap::elements_in_chunk) {}


ElementCacheMap::~ElementCacheMap() {
    if (element_eval_points_map_!=nullptr) {
        delete[] element_eval_points_map_;
    }
}


void ElementCacheMap::init(std::shared_ptr<EvalPoints> eval_points) {
    this->eval_points_ = eval_points;
    unsigned int ep_data_size = std::max(1.1, (double)eval_points_->max_size()) * CacheMapElementNumber::get();
    if (eval_point_data_.reserved_size() < ep_data_size) eval_point_data_.resize(ep_data_size);
    element_eval_points_map_ = new int [ep_data_size];
    for (unsigned int i=0; i<ep_data_size; ++i)
    	element_eval_points_map_[i] = ElementCacheMap::unused_point;
}


void ElementCacheMap::create_patch() {
    RevertableList<EvalPointData> eval_point_data_tmp = eval_point_data_;
    std::sort(eval_point_data_tmp.begin(), eval_point_data_tmp.end());
    eval_point_data_.reset();

    unsigned int last_region_idx = -1;
    unsigned int last_element_idx = -1;
    unsigned int i_pos=0; // position in eval_point_data_
    bool is_new_reg, is_new_elm;

    // Erase element data of previous step
    regions_starts_.reset();
    element_starts_.reset();
    element_to_map_.clear();
    element_to_map_bdr_.clear();
    std::fill(elm_idx_.begin(), elm_idx_.end(), ElementCacheMap::undef_elem_idx);

    for (auto it=eval_point_data_tmp.begin(); it!=eval_point_data_tmp.end(); ++it) {
        is_new_reg = (it->i_reg_ != last_region_idx);
        is_new_elm = is_new_reg || (it->i_element_ != last_element_idx);
        if (is_new_elm) {
            if (is_new_reg) {
                unsigned int last_eval_point = i_pos-1; // set size of block by SIMD size
                while (i_pos % ElementCacheMap::simd_size_double > 0) {
                	eval_point_data_.emplace_back( eval_point_data_[last_eval_point] );
                    i_pos++;
                }

                regions_starts_.emplace_back( element_starts_.temporary_size() );
                last_region_idx = it->i_reg_;
            }
			elm_idx_[element_starts_.temporary_size()] = it->i_element_;
            if (it->i_reg_ % 2 == 1) // bulk region > to element_to_map_ (bulk)
			    element_to_map_[it->i_element_] = element_starts_.temporary_size();
            else // boundary region to element_to_map_bdr_ (boundary)
                element_to_map_bdr_[it->i_element_] = element_starts_.temporary_size();
            element_starts_.emplace_back(i_pos);
            last_element_idx = it->i_element_;
        }
        eval_point_data_.emplace_back( *it );
        set_element_eval_point(element_starts_.temporary_size()-1, it->i_eval_point_, i_pos);
        i_pos++;
    }
    unsigned int last_eval_point = i_pos-1; // set size of block of last region by SIMD size
    while (i_pos % ElementCacheMap::simd_size_double > 0) {
        eval_point_data_.emplace_back( eval_point_data_[last_eval_point] );
        i_pos++;
    }

    regions_starts_.emplace_back( element_starts_.temporary_size() );
    element_starts_.emplace_back(i_pos);
    regions_starts_.make_permanent();
    element_starts_.make_permanent();
    eval_point_data_.make_permanent();
}


void ElementCacheMap::start_elements_update() {
	ready_to_reading_ = false;
}

void ElementCacheMap::finish_elements_update() {
	ready_to_reading_ = true;
}

