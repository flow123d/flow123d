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
 * @file    field_value_cache.impl.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_VALUE_CACHE_IMPL_HH_
#define FIELD_VALUE_CACHE_IMPL_HH_

#include "fields/field_value_cache.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/accessors.hh"

template<class elm_type, class return_type>
template<uint nRows, uint nCols>
Armor::Mat<elm_type, nRows, nCols> FieldValueCache<elm_type, return_type>::get_value(DHCellAccessor dh_cell,
        unsigned int subset_idx, unsigned int eval_points_idx) {

    ASSERT(dh_cell.element_cache_index() != ElementCacheMap::undef_elem_idx)(dh_cell.elm_idx());
    unsigned int points_per_element = this->subset_size(subset_idx) / n_cache_elements_;
    unsigned int subset_point_idx = eval_points_idx - eval_points_->subset_begin(subset_idx);
    return data_.get<nRows, nCols>(this->subset_begin(subset_idx) + dh_cell.element_cache_index() * points_per_element + subset_point_idx);
}



#endif /* FIELD_VALUE_CACHE_IMPL_HH_ */
