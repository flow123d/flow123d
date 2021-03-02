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

template<class Value>
typename Value::return_type ElementCacheMap::get_value(const FieldValueCache<typename Value::element_type> &field_cache,
        unsigned int elem_patch_idx, unsigned int eval_points_idx) const {
    ASSERT_EQ_DBG(Value::NRows_, field_cache.n_rows());
    ASSERT_EQ_DBG(Value::NCols_, field_cache.n_cols());
    unsigned int value_cache_idx = this->element_eval_point(elem_patch_idx, eval_points_idx);
    ASSERT_DBG(value_cache_idx != ElementCacheMap::undef_elem_idx);
    return Value::get_from_array(field_cache, value_cache_idx);
}



#endif /* FIELD_VALUE_CACHE_IMPL_HH_ */
