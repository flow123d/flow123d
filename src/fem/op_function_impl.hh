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
 * @file    op_function.hh
 * @brief   Store finite element reinit functions.
 * @author  David Flanderka
 */

#ifndef OP_FUNCTION_IMPL_HH_
#define OP_FUNCTION_IMPL_HH_

#include <Eigen/Dense>

#include "fem/patch_fe_values.hh"


template<unsigned int spacedim>
Scalar PatchOp<spacedim>::scalar_elem_value(uint point_idx) const {
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    return result_(0)( ppv.int_table_(1)(ppv.points_map_[point_idx]) );
}

template<unsigned int spacedim>
Vector PatchOp<spacedim>::vector_elem_value(uint point_idx) const {
    Vector val;
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    uint op_matrix_idx = ppv.int_table_(1)(ppv.points_map_[point_idx]);
    for (uint i=0; i<3; ++i)
        val(i) = result_(i)(op_matrix_idx);
    return val;
}

template<unsigned int spacedim>
Tensor PatchOp<spacedim>::tensor_elem_value(uint point_idx) const {
    Tensor val;
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    uint op_matrix_idx = ppv.int_table_(1)(ppv.points_map_[point_idx]);
    for (uint i=0; i<3; ++i)
        for (uint j=0; j<3; ++j)
            val(i,j) = result_(i+j*spacedim)(op_matrix_idx);
    return val;
}

template<unsigned int spacedim>
Scalar PatchOp<spacedim>::scalar_value(uint point_idx, uint i_dof) const {
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    return result_(i_dof)(ppv.points_map_[point_idx]);
}

template<unsigned int spacedim>
Vector PatchOp<spacedim>::vector_value(uint point_idx, uint i_dof) const {
    Vector val;
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<3; ++i)
        val(i) = result_(i + 3*i_dof)(op_matrix_idx);
    return val;
}

template<unsigned int spacedim>
Tensor PatchOp<spacedim>::tensor_value(uint point_idx, uint i_dof) const {
    Tensor val;
    PatchPointValues<spacedim> &ppv = patch_fe_->patch_point_vals_[bulk_side_][dim_-1];
    uint op_matrix_idx = ppv.points_map_[point_idx];
    for (uint i=0; i<9; ++i)
        val(i) = result_(i+9*i_dof)(op_matrix_idx);
    return val;
}


namespace Op {

namespace Bulk {

namespace El {

template<unsigned int spacedim>
void OpCoords<spacedim>::eval() {
    PatchPointValues<spacedim> &ppv = this->patch_fe_->patch_point_vals_[0][this->dim_-1];
    auto result = this->result_matrix();

    for (uint i_elm=0; i_elm<ppv.elem_list_.size(); ++i_elm)
        for (uint i_col=0; i_col<this->dim_+1; ++i_col)
            for (uint i_row=0; i_row<spacedim; ++i_row) {
                result(i_row, i_col)(i_elm) = ( *ppv.elem_list_[i_elm].node(i_col) )(i_row);
            }
}

template<unsigned int spacedim>
OpJac<spacedim>::OpJac(uint dim, PatchFEValues<spacedim> &pfev)
: PatchOp<spacedim>(dim, pfev, {spacedim, dim}, OpSizeType::elemOp)
{
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.template get< OpCoords<spacedim> >(dim) );
}

template<unsigned int op_dim, unsigned int spacedim>
OpInvJac<op_dim, spacedim>::OpInvJac(uint dim, PatchFEValues<3> &pfev)
: PatchOp<spacedim>(dim, pfev, {dim, spacedim}, OpSizeType::elemOp)
{
    ASSERT_EQ(this->dim_, op_dim);
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.template get< OpJac<spacedim> >(op_dim) );
}

template<unsigned int op_dim, unsigned int spacedim>
OpJacDet<op_dim, spacedim>::OpJacDet(uint dim, PatchFEValues<spacedim> &pfev)
: PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::elemOp)
{
    ASSERT_EQ(this->dim_, dim);
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.template get< OpJac<spacedim> >(dim) );
}

} // end of namespace Op::Bulk::El

namespace Pt {

template<unsigned int spacedim>
OpWeights<spacedim>::OpWeights(uint dim, PatchFEValues<spacedim> &pfev)
: PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::fixedSizeOp)
{
    this->bulk_side_ = 0;

    // create result vector of weights operation in assembly arena
    const std::vector<double> &point_weights_vec = pfev.get_bulk_quadrature(dim)->get_weights();
    this->allocate_result(point_weights_vec.size(), pfev.asm_arena());
    for (uint i=0; i<point_weights_vec.size(); ++i)
        this->result_(0)(i) = point_weights_vec[i];
}

template<unsigned int op_dim, unsigned int spacedim>
OpJxW<op_dim, spacedim>::OpJxW(uint dim, PatchFEValues<spacedim> &pfev)
: PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::pointOp)
{
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.template get< OpWeights<spacedim> >(dim) );
    this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpJacDet<op_dim, spacedim> >(dim) );
}

template<unsigned int op_dim, unsigned int spacedim>
OpRefGradScalar<op_dim, spacedim>::OpRefGradScalar(uint dim, PatchFEValues<spacedim> &pfev, uint component_idx)
: PatchOp<spacedim>(dim, pfev, {dim, 1}, OpSizeType::fixedSizeOp)
{
    this->bulk_side_ = 0;
    auto fe_component = pfev.template fe_comp<op_dim>(component_idx);
    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

    uint n_points = pfev.get_bulk_quadrature(op_dim)->size();
    uint n_dofs = fe_component->n_dofs();
    this->n_dofs_ = n_dofs;
    this->shape_[1] = n_dofs;

    std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<op_dim>(pfev.get_bulk_quadrature(op_dim), fe_component);
    this->allocate_result(n_points, pfev.asm_arena());
    auto ref_scalar_value = this->result_matrix();
    for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row) {
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                ref_scalar_value(i_row, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i_row);
    }
}

template<unsigned int op_dim, unsigned int spacedim>
OpGradScalarShape<op_dim, spacedim>::OpGradScalarShape(uint dim, PatchFEValues<spacedim> &pfev, uint component_idx)
: PatchOp<spacedim>(dim, pfev, {spacedim, 1}, OpSizeType::pointOp)
{
    this->bulk_side_ = 0;
    auto fe_component = pfev.template fe_comp<op_dim>(component_idx);
    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

    uint n_dofs = fe_component->n_dofs();
    this->n_dofs_ = n_dofs;
    this->shape_[1] = n_dofs;

    this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpInvJac<op_dim, spacedim> >(dim) );
    this->input_ops_.push_back( pfev.template get< OpRefGradScalar<op_dim, spacedim> >(dim, component_idx) );
}

} // end of namespace Op::Bulk::Pt

} // end of namespace Op::Bulk

namespace Side {

namespace El {

} // end of namespace Op::Side::El

namespace Pt {

} // end of namespace Op::Side::Pt

} // end of namespace Side

} // end of namespace Op


// explicit instantiation of template classes
template class Op::Bulk::El::OpCoords<3>;
template class Op::Bulk::El::OpJac<3>;
template class Op::Bulk::El::OpInvJac<1, 3>;
template class Op::Bulk::El::OpInvJac<2, 3>;
template class Op::Bulk::El::OpInvJac<3, 3>;
template class Op::Bulk::El::OpJacDet<1, 3>;
template class Op::Bulk::El::OpJacDet<2, 3>;
template class Op::Bulk::El::OpJacDet<3, 3>;
template class Op::Bulk::Pt::OpJxW<1, 3>;
template class Op::Bulk::Pt::OpJxW<2, 3>;
template class Op::Bulk::Pt::OpJxW<3, 3>;
template class Op::Bulk::Pt::OpRefGradScalar<1, 3>;
template class Op::Bulk::Pt::OpRefGradScalar<2, 3>;
template class Op::Bulk::Pt::OpRefGradScalar<3, 3>;
template class Op::Bulk::Pt::OpGradScalarShape<1, 3>;
template class Op::Bulk::Pt::OpGradScalarShape<2, 3>;
template class Op::Bulk::Pt::OpGradScalarShape<3, 3>;

template class PatchOp<3>;


#endif /* OP_FUNCTION_IMPL_HH_ */
