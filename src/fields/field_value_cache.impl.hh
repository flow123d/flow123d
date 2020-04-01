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

    ASSERT(dh_cell.element_cache_index() != ElementCacheMap::undef_elem_idx)(dh_cell.elm_idx());
    int value_cache_idx = map.get_field_value_cache_index(dh_cell.element_cache_index(), eval_points_idx);
    ASSERT_GE(value_cache_idx, 0);
    return Value::get_from_array(data_, value_cache_idx);
}



#endif /* FIELD_VALUE_CACHE_IMPL_HH_ */
