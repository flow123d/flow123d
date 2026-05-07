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
 * @file    op_accessors.hh
 * @brief   Declares accessors to FE operations.
 * @author  David Flanderka
 */

#ifndef OP_ACCESSORS_IMPL_HH_
#define OP_ACCESSORS_IMPL_HH_

#include "fem/patch_fe_values.hh"
#include "fem/op_accessors.hh"
#include "fem/patch_op.hh"


template <class ValueType>
inline ValueType ElQ<ValueType>::operator()(const BulkPoint &point) const {
    return patch_op_->elem_value<ValueType>( point.value_cache_idx() );
}

template <class ValueType>
inline ValueType ElQ<ValueType>::operator()(const SidePoint &point) const {
    return patch_op_->elem_value<ValueType>( point.value_cache_idx() );
}

template <class ValueType>
inline ValueType FeQ<ValueType>::operator()(const BulkPoint &point) const {
    ASSERT_PTR(patch_op_bulk_);
    return patch_op_bulk_->point_value<ValueType>(point.value_cache_idx(), i_shape_fn_idx_);
}

template <class ValueType>
inline ValueType FeQ<ValueType>::operator()(const SidePoint &point) const {
    ASSERT_PTR(patch_op_side_);
    return patch_op_side_->point_value<ValueType>(point.value_cache_idx(), i_shape_fn_idx_);
}

#endif /* OP_ACCESSORS_IMPL_HH_ */
