/*!
 *
ÄŹÂ»ĹĽ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
#include "fem/patch_point_values.hh"
#include "fem/patch_op.hh"
#include "quadrature/quadrature.hh"
#include <cmath>
#include <initializer_list>

namespace Op {

/**
 * Holds common functionality of patch operations.
 */
template<unsigned int spacedim = 3>
class FuncHelper {
public:
    /**
     * Copy reduced data from 'source' to 'target' ArenaVec. Mapping of reduced data is given by 'ppv' data.
     */
    static void fill_reduce_element_data_vec(PatchPointValues<spacedim> &ppv, ArenaVec<double> &source, ArenaVec<double> &target) {
        for (uint i_el=0; i_el<ppv.n_mesh_items(); ++i_el) {
            target( i_el ) = source( ppv.int_table_(patch_elem_on_domain)(i_el) );
        }
    }

    /**
     * Copy ref_side_on_quads data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_quads_data(PatchOp<spacedim> &op, unsigned int i_comp,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &target) {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_dofs = op.n_dofs();
        uint n_sides = ppv.n_mesh_items();               // number of sides on patch
        uint n_patch_points = n_sides * op.quad_size();  // number of points on patch
        uint n_comp = op.shape()[0];                     // number of components of result

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                target(i_comp,i_dof)(i_pt) = source(ppv.int_table_(ref_side_on_quads)(i_pt),n_comp*i_dof+i_comp)(i_pt / n_sides);
        }
    }

    /**
     * Copy ref_side_on_sides vector data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_sides_vector_data(PatchOp<spacedim> &op, unsigned int dim,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> &target)
    {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_sides = ppv.n_mesh_items();
        for (uint i_dof=0; i_dof<op.n_dofs(); ++i_dof) {
            for (uint i_pt=0; i_pt<op.quad_size(); ++i_pt) {
                uint i_begin = i_pt * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<dim; ++i_c) {
                        target(i_c, i_dof)(i_begin + i_sd) = source(ppv.int_table_(ref_side_on_sides)(i_sd), i_dof*dim+i_c)(i_pt);
                    }
                }
            }
        }
    }

    /**
     * Copy ref_side_on_sides tensor data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_sides_tensor_data(PatchOp<spacedim> &op, unsigned int i_dof, unsigned int shape_m, unsigned int shape_n,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
			Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &target)
    {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_sides = ppv.n_mesh_items();
        uint n_points = source(0).data_size(); // number of ref points

        for (uint i_pt=0; i_pt<n_points; ++i_pt) {
            uint i_begin = i_pt * n_sides;
            for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                for (uint i_dim=0; i_dim<shape_m; ++i_dim) {
                    for (uint i_c=0; i_c<shape_n; ++i_c) {
                        target(i_dim, i_c)(i_begin + i_sd) = source(ppv.int_table_(ref_side_on_sides)(i_sd)*shape_m+i_dim, shape_n*i_dof+i_c)(i_pt);
                    }
                }
            }
        }
    }

    /**
     * Copy patch_elem_on_domain data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_patch_elem_on_domain_data(PatchOp<spacedim> &op, unsigned int dim,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &target) {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_elems = ppv.n_mesh_items();          // number of elements on patch
        uint n_elem_points = op.quad_size();        // number of points on element

        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            target(i_c) = ArenaVec<double>( n_elems*n_elem_points, op.patch_arena() );
            for (uint i_el=0; i_el<n_elems; ++i_el)
                for (uint i_pt=0; i_pt<n_elem_points; ++i_pt) {
                    target(i_c)( i_pt*n_elems + i_el ) = source(i_c)( ppv.int_table_(patch_elem_on_domain)(i_el) );
                }
        }
    }
private:
    /// Forbidden constructor
    FuncHelper() {}
};


/**
 * @brief Common types for accumulating finite element shape function values.
 */
template <typename Derived>
class ShapeFnMatrixBase {
protected:
    using Matrix = Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixMap = Eigen::Map<Matrix>;
    using DofValues = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>;

public:
    /**
     * @brief Initialize result matrix to zero values.
     */
    static void initialize_reference(MatrixMap result) {
        Derived::initialize_reference(result);
    }

    /**
     * @brief Store value of shape function to result matrix.
     *
     * Position in result matrix is given by index of side, row, column, DOF and quadrature point.
     */
    template<unsigned int N>
    static void store_reference(MatrixMap result, unsigned int side, unsigned int row,
                                unsigned int col, unsigned int dof, unsigned int point, double value)
    {
        Derived::template store_reference<N>(result, side, row, col, dof, point, value);
    }

    /**
     * @brief Expand reference values of shape function on bulk domain ​​for all quadrature points on patch
     */
    static void expand_bulk(PatchOp<3> &op) {
        Derived::expand_bulk(op);
    }

    /**
     * @brief Expand reference values of shape function on side domain ​​for all quadrature points on
     */
    static void expand_side(PatchOp<3> &op) {
        Derived::expand_side(op);
    }

    /**
     * @brief Accumulate result over all shape functions (DOFs).
     */
    template<class Domain>
    static void integrate_dofs(MatrixMap result, PatchOp<3> &shape_op, const DofValues &data, Quadrature &quad)
    {
        Derived::template integrate_dofs<Domain>(result, shape_op, data, quad);
    }
};

/// Accumulate all components of each shape function matrix.
class ShapeFnMatrix : public ShapeFnMatrixBase<ShapeFnMatrix> {
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

    static void expand_bulk(PatchOp<3> &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        auto result = op.result_matrix();
        ArenaVec<double> elements(op.ppv().n_mesh_items(), op.patch_arena());
        for (unsigned int i=0; i<elements.data_size(); ++i) elements(i) = 1.0;
        ArenaOVec<double> element_values(elements);
        uint n_rows = op.shape()[0], n_cols = op.shape()[1];

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_values(n_rows, n_cols);
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof) {
            for (unsigned int c=0; c<n_rows*n_cols; ++c)
                ref_values(c) = ArenaOVec<double>(source(dof*n_rows*n_cols+c));
            Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> values = element_values * ref_values;
            for (unsigned int c=0; c<n_rows*n_cols; ++c)
                result(dof*n_rows*n_cols+c) = values(c).get_vec();
        }
    }

    static void expand_side(PatchOp<3> &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        for (unsigned int dof=0; dof<op.n_dofs(); ++dof) {
            auto result = op.result_sub_matrix(dof);
            Op::FuncHelper<3>::copy_ref_side_on_sides_tensor_data(op, dof, op.shape()[0], op.shape()[1], source, result);
        }
    }

    template<class Domain>
    static void integrate_dofs(MatrixMap result, PatchOp<3> &shape_op, const DofValues &data, Quadrature &)
    {
        for (unsigned int i_dof=0; i_dof<shape_op.n_dofs(); ++i_dof) {
            MatrixMap input_submat = shape_op.result_sub_matrix(i_dof);
            result = result + data(i_dof) * input_submat;
        }
    }
};

/// Accumulate scalar shape values into their respective result components.
class ShapeFnCompressed : public ShapeFnMatrixBase<ShapeFnCompressed> {
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

    static void expand_bulk(PatchOp<3> &op)
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

    static void expand_side(PatchOp<3> &op)
    {
        auto source = op.input_ops(0)->result_matrix();
        auto result = op.result_matrix();
        Op::FuncHelper<3>::copy_ref_side_on_quads_data(op, 0, source, result);
    }

    template<class Domain>
    static void integrate_dofs(MatrixMap result, PatchOp<3> &shape_op, const DofValues &data, Quadrature &quad)
    {
        unsigned int n_bulk_points = Domain::n_bulk_points(quad);
        for (unsigned int i_dof=0; i_dof<shape_op.n_dofs(); ++i_dof) {
            result(i_dof/n_bulk_points) = result(i_dof/n_bulk_points)
                + data(i_dof) * shape_op.result_sub_matrix(i_dof)(0);
        }
    }
};


}  // end of namespace Op

#endif /* FEM_SHAPE_FN_MATRIX_HH_ */
