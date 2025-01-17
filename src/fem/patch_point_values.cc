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
 * @file    patch_point_values.cc
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#include "fem/patch_point_values.hh"


template<unsigned int spacedim>
void PatchPointValues<spacedim>::create_zero_operations(std::vector<PatchOp<spacedim> *> &ref_ops) {
    operations_.resize(ref_ops.size(), nullptr);
    for (uint i_op = 0; i_op < ref_ops.size(); ++i_op ) {
        auto *op = ref_ops[i_op];
        if (op == nullptr) continue;

        auto *new_op = make_fe_op(i_op, {op->shape()[0], op->shape()[1]}, &common_reinit::op_base, op->n_dofs(), op->size_type());
        new_op->allocate_const_result(patch_fe_data_.zero_vec_);
    }
}


template class PatchPointValues<3>;
template class PatchOp<3>;
