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

#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/eval_points.hh"


template<class Value>
const unsigned int FieldValueCache<Value>::n_cached_elements = 20;

template<class Value>
FieldValueCache<Value>::FieldValueCache()
: data_(nullptr) {}

template<class Value>
FieldValueCache<Value>::FieldValueCache(EvalPoints eval_points) {
	unsigned int size = n_cached_elements * eval_points.size() * Value::NRows_ * Value::NCols_;
	data_ = new double[size];
}

template<class Value>
FieldValueCache<Value>::~FieldValueCache() {
	if (data_!=nullptr) delete [] data_;
}

template<class Value>
void FieldValueCache<Value>::mark_used(EvalSubset sub_quad) {
    auto local_point_vec = sub_quad.get_point_indices(0);
    for (auto p_idx : local_point_vec) used_points_.insert(p_idx);
}


/******************************************************************************
 * Explicit instantiation of templates
 */

//template class FieldValueCache<FieldValue<0>::Enum >;
//template class FieldValueCache<FieldValue<0>::Integer >;
template class FieldValueCache<FieldValue<0>::Scalar >;
template class FieldValueCache<FieldValue<3>::VectorFixed >;
template class FieldValueCache<FieldValue<3>::TensorFixed >;
