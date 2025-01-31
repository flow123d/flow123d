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

#include "fem/patch_point_values.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_accessors.hh"
#include "fem/op_function.hh"


template <class ValueType>
ValueType ElQ<ValueType>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->scalar_elem_value(value_cache_idx);
//    return patch_point_vals_->scalar_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->vector_elem_value(value_cache_idx);
//    return patch_point_vals_->vector_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->tensor_elem_value(value_cache_idx);
//    return patch_point_vals_->tensor_elem_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType ElQ<ValueType>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->scalar_elem_value(value_cache_idx);
//    return patch_point_vals_->scalar_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->vector_elem_value(value_cache_idx);
//    return patch_point_vals_->vector_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_->tensor_elem_value(value_cache_idx);
//    return patch_point_vals_->tensor_elem_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(const BulkPoint &point) const {
    ASSERT_PTR(patch_op_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_bulk_->scalar_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_bulk_->scalar_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Vector FeQ<Vector>::operator()(const BulkPoint &point) const {
    ASSERT_PTR(patch_op_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_bulk_->vector_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_bulk_->vector_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Tensor FeQ<Tensor>::operator()(const BulkPoint &point) const {
    ASSERT_PTR(patch_op_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_bulk_->tensor_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_bulk_->tensor_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(const SidePoint &point) const {
    ASSERT_PTR(patch_op_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_side_->scalar_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_side_->scalar_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Vector FeQ<Vector>::operator()(const SidePoint &point) const {
    ASSERT_PTR(patch_op_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_side_->vector_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_side_->vector_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Tensor FeQ<Tensor>::operator()(const SidePoint &point) const {
    ASSERT_PTR(patch_op_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_op_side_->tensor_value(value_cache_idx, i_shape_fn_idx_);
//    return patch_point_vals_side_->tensor_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}


#endif /* OP_ACCESSORS_IMPL_HH_ */
