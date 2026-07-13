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


template<unsigned int spacedim>
template <unsigned int elemdim>
inline arma::mat::fixed<spacedim, elemdim> PatchOp<spacedim>::elem_matrix_value(uint point_idx) const {
	arma::mat::fixed<spacedim, elemdim> val;
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    uint op_matrix_idx = ppv.int_table_(domain_on_quads)(ppv.points_map_[point_idx]);
    for (uint i=0; i<spacedim; ++i)
        for (uint j=0; j<elemdim; ++j)
            val(i,j) = result_(i+j*spacedim)(op_matrix_idx);
    return val;
}

template<>
template<>
inline Scalar PatchOp<3>::elem_value<Scalar>(uint point_idx) const {
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    return result_(0)( ppv.int_table_(domain_on_quads)(ppv.points_map_[point_idx]) );
}

template<>
template<>
inline Vector PatchOp<3>::elem_value<Vector>(uint point_idx) const {
    Vector val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    uint op_matrix_idx = ppv.int_table_(domain_on_quads)(ppv.points_map_[point_idx]);
    for (uint i=0; i<3; ++i)
        val(i) = result_(i)(op_matrix_idx);
    return val;
}

template<>
template<>
inline arma::mat::fixed<3,0> PatchOp<3>::elem_value<arma::mat::fixed<3,0>>(uint point_idx) const {
    return this->elem_matrix_value<0>(point_idx);
}

template<>
template<>
inline arma::mat::fixed<3,1> PatchOp<3>::elem_value<arma::mat::fixed<3,1>>(uint point_idx) const {
    return this->elem_matrix_value<1>(point_idx);
}

template<>
template<>
inline arma::mat::fixed<3,2> PatchOp<3>::elem_value<arma::mat::fixed<3,2>>(uint point_idx) const {
    return this->elem_matrix_value<2>(point_idx);
}

template<>
template<>
inline arma::mat::fixed<3,3> PatchOp<3>::elem_value<arma::mat::fixed<3,3>>(uint point_idx) const {
//    arma::mat::fixed<3,3> val;
//    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
//    uint op_matrix_idx = ppv.int_table_(domain_on_quads)(ppv.points_map_[point_idx]);
//    for (uint i=0; i<3; ++i)
//        for (uint j=0; j<3; ++j)
//            val(i,j) = result_(i+j*3)(op_matrix_idx);
//    return val;
    return this->elem_matrix_value<3>(point_idx);
}

template<>
template<>
inline Scalar PatchOp<3>::point_value<Scalar>(uint point_idx, uint i_dof) const {
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    return result_(i_dof)(ppv.points_map_[point_idx]);
}

template<>
template<>
inline Vector PatchOp<3>::point_value<Vector>(uint point_idx, uint i_dof) const {
    Vector val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<3; ++i)
        val(i) = result_(i + 3*i_dof)(op_matrix_idx);
    return val;
}

template<>
template<>
inline Tensor PatchOp<3>::point_value<Tensor>(uint point_idx, uint i_dof) const {
    Tensor val;
    PatchPointValues<3> &ppv = patch_fe_->patch_point_vals_[domain_][dim_];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<9; ++i)
        val(i) = result_(i+9*i_dof)(op_matrix_idx);
    return val;
}

template <>
template <>
inline unsigned int PatchOp<3>::point_value<unsigned int>(FMT_UNUSED uint point_idx, FMT_UNUSED uint i_dof) const {
    return 0;
}

template <>
template <>
inline int PatchOp<3>::point_value<int>(FMT_UNUSED uint point_idx, FMT_UNUSED uint i_dof) const {
    return 0;
}




// explicit instantiation of template classes

template class PatchOp<3>;


#endif /* OP_FUNCTION_IMPL_HH_ */
