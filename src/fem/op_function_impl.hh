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

#include "fem/op_function.hh"
#include "fem/patch_fe_values.hh"

namespace Op {

namespace Bulk {

namespace El {

OpJac::OpJac(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, {3, dim}, OpSizeType::elemOp)
{
    this->input_ops_.push_back( pfev.get< OpCoords >(dim) );
}

OpInvJac::OpInvJac(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, {dim, 3}, OpSizeType::elemOp)
{
    this->input_ops_.push_back( pfev.get< OpJac >(dim) );
}

OpJacDet::OpJacDet(uint dim, PatchFEValues<3> &pfev)
: PatchOp<3>(dim, {1}, OpSizeType::elemOp)
{
    this->input_ops_.push_back( pfev.get< OpJac >(dim) );
}

} // end of namespace Op::Bulk::El

namespace Pt {

//OpRefGradScalar::OpRefGradScalar(uint dim, PatchFEValues<3> &pfev, uint n_dofs)
//: PatchOp<3>(dim, {dim, n_dofs}, OpSizeType::fixedOp, n_dofs)
//{
//    uint n_points = patch_point_vals_->get_quadrature()->size(); // get quadrature size
//
//    std::vector<std::vector<arma::mat> > ref_shape_grads = this->ref_shape_gradients_bulk(patch_point_vals_->get_quadrature(), fe_component); //move to PatchOp
//    this->allocate_result(n_points, pfev.asm_arena());
//    auto ref_scalar_value = this->result_matrix();
//    for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row) {
//        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
//            for (uint i_p=0; i_p<n_points; ++i_p)
//                ref_scalar_value(i_row, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i_row);
//    }
//
//}

} // end of namespace Op::Bulk::Pt

} // end of namespace Op::Bulk

namespace Side {

namespace El {

} // end of namespace Op::Side::El

namespace Pt {

} // end of namespace Op::Side::Pt

} // end of namespace Side

} // end of namespace Op


#endif /* OP_FUNCTION_IMPL_HH_ */
