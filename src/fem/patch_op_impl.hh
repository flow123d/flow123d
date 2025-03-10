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
 * @file    patch_op_impl.hh
 * @brief   Store finite element reinit functions.
 * @author  David Flanderka
 */

#ifndef OP_FUNCTION_IMPL_HH_
#define OP_FUNCTION_IMPL_HH_

#include <Eigen/Dense>

#include "fem/patch_op.hh"
#include "fem/patch_fe_values.hh"


template<>
template<>
Scalar PatchOp<3>::elem_value<Scalar>(uint point_idx) const {
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    return result_(0)( ppv.int_table_(1)(ppv.points_map_[point_idx]) );
}

template<>
template<>
Vector PatchOp<3>::elem_value<Vector>(uint point_idx) const {
    Vector val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    uint op_matrix_idx = ppv.int_table_(1)(ppv.points_map_[point_idx]);
    for (uint i=0; i<3; ++i)
        val(i) = result_(i)(op_matrix_idx);
    return val;
}

template<>
template<>
Tensor PatchOp<3>::elem_value<Tensor>(uint point_idx) const {
    Tensor val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    uint op_matrix_idx = ppv.int_table_(1)(ppv.points_map_[point_idx]);
    for (uint i=0; i<3; ++i)
        for (uint j=0; j<3; ++j)
            val(i,j) = result_(i+j*3)(op_matrix_idx);
    return val;
}

template<>
template<>
Scalar PatchOp<3>::point_value<Scalar>(uint point_idx, uint i_dof) const {
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    return result_(i_dof)(ppv.points_map_[point_idx]);
}

template<>
template<>
Vector PatchOp<3>::point_value<Vector>(uint point_idx, uint i_dof) const {
    Vector val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<3; ++i)
        val(i) = result_(i + 3*i_dof)(op_matrix_idx);
    return val;
}

template<>
template<>
Tensor PatchOp<3>::point_value<Tensor>(uint point_idx, uint i_dof) const {
    Tensor val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_-1];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<9; ++i)
        val(i) = result_(i+9*i_dof)(op_matrix_idx);
    return val;
}




// explicit instantiation of template classes

template class PatchOp<3>;


#endif /* OP_FUNCTION_IMPL_HH_ */
