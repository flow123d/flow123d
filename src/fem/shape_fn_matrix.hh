/*!
 *
ď»ż * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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

#ifndef FEM_SHAPE_FN_MATRIX_HH_
#define FEM_SHAPE_FN_MATRIX_HH_

#include <Eigen/Dense>
#include "fem/arena_vec.hh"
#include <cmath>
#include <initializer_list>

namespace Op {
template<unsigned int spacedim> class FuncHelper;
}


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
    /// Storage dimensions, including reference values for all sides if requested.
    template<unsigned int N, unsigned int NSides = 1>
    struct StorageShape {
        static constexpr std::initializer_list<unsigned int> value = {NSides*N, N};
    };

    static void initialize_reference(MatrixMap) {}

    template<unsigned int N>
    static void store_reference(MatrixMap result, unsigned int side, unsigned int row,
                                unsigned int col, unsigned int dof, unsigned int point, double value)
    {
        result(side*N+row, N*dof+col)(point) = value;
    }

    template<unsigned int N, class ShapeOp>
    static void expand_bulk(ShapeOp &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        auto result = op.result_matrix();
        ArenaVec<double> elements(op.ppv().n_mesh_items(), op.patch_arena());
        for (unsigned int i=0; i<elements.data_size(); ++i) elements(i) = 1.0;
        ArenaOVec<double> element_values(elements);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_values(N, N);
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof) {
            for (unsigned int c=0; c<N*N; ++c)
                ref_values(c) = ArenaOVec<double>(source(dof*N*N+c));
            Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> values = element_values * ref_values;
            for (unsigned int c=0; c<N*N; ++c)
                result(dof*N*N+c) = values(c).get_vec();
        }
    }

    template<unsigned int spacedim, unsigned int N, class ShapeOp>
    static void expand_side(ShapeOp &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof) {
            auto result = op.result_sub_matrix(dof);
            Op::FuncHelper<spacedim>::copy_ref_side_on_sides_tensor_data(op, dof, N, N, source, result);
        }
    }

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
    /// One scalar per DOF, with one row per reference side.
    template<unsigned int N, unsigned int NSides = 1>
    struct StorageShape {
        static constexpr std::initializer_list<unsigned int> value = {NSides};
    };

    // A shape function may vanish at a quadrature point. Such values must remain zero.
    static void initialize_reference(MatrixMap result)
    {
        for (Eigen::Index c=0; c<result.size(); ++c)
            for (unsigned int point=0; point<result(c).data_size(); ++point)
                result(c)(point) = 0.0;
    }

    template<unsigned int N>
    static void store_reference(MatrixMap result, unsigned int side, unsigned int,
                                unsigned int, unsigned int dof, unsigned int point, double value)
    {
        if (std::fabs(value) > 1e-13) result(side, dof)(point) = value;
    }

    template<unsigned int N, class ShapeOp>
    static void expand_bulk(ShapeOp &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        auto result = op.result_matrix();
        ArenaVec<double> elements(op.ppv().n_mesh_items(), op.patch_arena());
        for (unsigned int i=0; i<elements.data_size(); ++i) elements(i) = 1.0;
        ArenaOVec<double> element_values(elements);

        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> ref_values(op.n_dofs());
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof)
            ref_values(dof) = ArenaOVec<double>(source(dof));
        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> values = element_values * ref_values;
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof)
            result(dof) = values(dof).get_vec();
    }

    template<unsigned int spacedim, unsigned int N, class ShapeOp>
    static void expand_side(ShapeOp &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        auto result = op.result_matrix();
        Op::FuncHelper<spacedim>::copy_ref_side_on_quads_data(op, 0, source, result);
    }

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

#endif /* FEM_SHAPE_FN_MATRIX_HH_ */
