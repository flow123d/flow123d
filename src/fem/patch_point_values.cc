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


// template specialization
template<>
void side_reinit::elop_sd_jac<1>(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
}

template<>
void side_reinit::elop_sd_jac_det<1>(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
    auto &op = operations[FeSide::SideOps::opSdJacDet];
    ArrayDbl &result_vec = op_results( op.result_row() );
    for (uint i=0;i<result_vec.size(); ++i) {
        result_vec(i) = 1.0;
    }
}


template class PatchPointValues<3>;
template class ElOp<3>;
template class FeBulk::PatchPointValues<3>;
template class FeSide::PatchPointValues<3>;
