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

template<class elm_type>
template<class Value>
typename Value::return_type FieldValueCache<elm_type>::get_value(const ElementCacheMap &map,
        const DHCellAccessor &dh_cell, unsigned int eval_points_idx) {
    static_assert( std::is_same<elm_type, typename Value::element_type>::value, "Wrong element type.");

    ASSERT_EQ_DBG(Value::NRows_, data_.n_rows());
    ASSERT_EQ_DBG(Value::NCols_, data_.n_cols());
    unsigned int value_cache_idx = map.element_eval_point(map.position_in_cache(dh_cell.elm().mesh_idx()), eval_points_idx);
    ASSERT_DBG(value_cache_idx != ElementCacheMap::undef_elem_idx);
    return Value::get_from_array(data_, value_cache_idx);
}


template<class elm_type>
template<class Value>
typename Value::return_type FieldValueCache<elm_type>::get_value(const ElementCacheMap &map,
        const ElementAccessor<3> elm, unsigned int eval_points_idx) {
    static_assert( std::is_same<elm_type, typename Value::element_type>::value, "Wrong element type.");

    ASSERT_EQ_DBG(Value::NRows_, data_.n_rows());
    ASSERT_EQ_DBG(Value::NCols_, data_.n_cols());
    unsigned int value_cache_idx = map.element_eval_point(map.position_in_cache(elm.mesh_idx()), eval_points_idx);
    ASSERT_DBG(value_cache_idx != ElementCacheMap::undef_elem_idx);
    return Value::get_from_array(data_, value_cache_idx);
}



#endif /* FIELD_VALUE_CACHE_IMPL_HH_ */
