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
 * @file    shape_fn_matrix.hh
 * @brief
 */

#ifndef SHAPE_FN_MATRIX_HH_
#define SHAPE_FN_MATRIX_HH_


#include <Eigen/Dense>
#include "fem/arena_vec.hh"

/// Common types for accumulating finite element shape function values.
class ShapeFnMatrixBase {
protected:
    using Matrix = Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixMap = Eigen::Map<Matrix>;
    using DofValues = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>;
};

/// Accumulate all components of each shape function matrix.
class ShapeFnMatrix : public ShapeFnMatrixBase {
public:
    template<class Domain, class ShapeOp, class Quad>
    static void eval(MatrixMap result, ShapeOp &shape_op, const DofValues &data, Quad &)
    {
        for (unsigned int i_dof=0; i_dof<shape_op.n_dofs(); ++i_dof) {
            MatrixMap input_submat = shape_op.result_sub_matrix(i_dof);
            result = result + data(i_dof) * input_submat;
        }
    }
};

/// Accumulate scalar shape values into their respective result components.
class ShapeFnCompressed : public ShapeFnMatrixBase {
public:
    template<class Domain, class ShapeOp, class Quad>
    static void eval(MatrixMap result, ShapeOp &shape_op, const DofValues &data, Quad &quad)
    {
        unsigned int n_bulk_points = Domain::n_bulk_points(quad);
        for (unsigned int i_dof=0; i_dof<shape_op.n_dofs(); ++i_dof) {
            result(i_dof/n_bulk_points) = result(i_dof/n_bulk_points)
                + data(i_dof) * shape_op.result_sub_matrix(i_dof)(0);
        }
    }
};

#endif /* SHAPE_FN_MATRIX_HH_ */
