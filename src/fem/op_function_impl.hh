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

OpCoords::OpCoords(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, pfev, {3, dim+1}, OpSizeType::elemOp)
{
    this->bulk_side_ = 0;
    pfev.patch_point_vals_[0][dim-1].op_el_coords_ = this;
}

OpJac::OpJac(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, pfev, {3, dim}, OpSizeType::elemOp)
{
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.get< OpCoords >(dim) );
}

template<unsigned int op_dim>
OpInvJac<op_dim>::OpInvJac(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, pfev, {dim, 3}, OpSizeType::elemOp)
{
    ASSERT_EQ(this->dim_, op_dim);
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.get< OpJac >(op_dim) );
}

template<unsigned int op_dim>
OpJacDet<op_dim>::OpJacDet(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, pfev, {1}, OpSizeType::elemOp)
{
    ASSERT_EQ(this->dim_, dim);
    this->bulk_side_ = 0;
    this->input_ops_.push_back( pfev.get< OpJac >(dim) );
}

} // end of namespace Op::Bulk::El

namespace Pt {

template<unsigned int op_dim>
OpRefGradScalar<op_dim>::OpRefGradScalar(uint dim, PatchFEValues<3> &pfev, uint component_idx)
: PatchOp<3>(dim, pfev, {dim, 1}, OpSizeType::fixedSizeOp)
{
    this->bulk_side_ = 0;
    auto fe_component = pfev.fe_comp<op_dim>(component_idx);
    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

    uint n_points = pfev.get_bulk_quadrature(op_dim)->size();
    uint n_dofs = fe_component->n_dofs();
    this->n_dofs_ = n_dofs;
    this->shape_[1] = n_dofs;

    std::vector<std::vector<arma::mat> > ref_shape_grads = this->ref_shape_gradients_bulk<op_dim>(pfev.get_bulk_quadrature(op_dim), fe_component);
    this->allocate_result(n_points, pfev.asm_arena());
    auto ref_scalar_value = this->result_matrix();
    for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row) {
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                ref_scalar_value(i_row, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i_row);
    }
}

template<unsigned int op_dim>
OpGradScalarShape<op_dim>::OpGradScalarShape(uint dim, PatchFEValues<3> &pfev, uint component_idx)
: PatchOp<3>(dim, pfev, {3, 1}, OpSizeType::pointOp)
{
    this->bulk_side_ = 0;
    auto fe_component = pfev.fe_comp<op_dim>(component_idx);
    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

    uint n_dofs = fe_component->n_dofs();
    this->n_dofs_ = n_dofs;
    this->shape_[1] = n_dofs;

    this->input_ops_.push_back( pfev.get< Op::Bulk::El::OpInvJac<op_dim> >(dim) );
    this->input_ops_.push_back( pfev.get< OpRefGradScalar<op_dim> >(dim, component_idx) );
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
template class Op::Bulk::El::OpInvJac<1>;
template class Op::Bulk::El::OpInvJac<2>;
template class Op::Bulk::El::OpInvJac<3>;
template class Op::Bulk::El::OpJacDet<1>;
template class Op::Bulk::El::OpJacDet<2>;
template class Op::Bulk::El::OpJacDet<3>;
template class Op::Bulk::Pt::OpRefGradScalar<1>;
template class Op::Bulk::Pt::OpRefGradScalar<2>;
template class Op::Bulk::Pt::OpRefGradScalar<3>;
template class Op::Bulk::Pt::OpGradScalarShape<1>;
template class Op::Bulk::Pt::OpGradScalarShape<2>;
template class Op::Bulk::Pt::OpGradScalarShape<3>;

template class PatchOp<3>;


#endif /* OP_FUNCTION_IMPL_HH_ */
