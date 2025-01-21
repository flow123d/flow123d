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

#ifndef OP_FUNCTION_HH_
#define OP_FUNCTION_HH_

#include <Eigen/Dense>

#include "fem/eigen_tools.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"
#include "mesh/ref_element.hh"



template<unsigned int spacedim> class PatchOp;


/// Type for conciseness
using ReinitFunction = std::function<void(PatchOp<3> *, IntTableArena &)>;


/// Distinguishes operations by type and size of output rows
enum OpSizeType
{
	elemOp,      ///< operation is evaluated on elements or sides
	pointOp,     ///< operation is evaluated on quadrature points
	fixedSizeOp  ///< operation has fixed size and it is filled during initialization
};



/**
 * @brief Class represents element or FE operations.
 */
template<unsigned int spacedim = 3>
class PatchOp {
public:
    /**
     * Constructor
     *
     * Set all data members.
     */
    PatchOp(uint dim, std::initializer_list<uint> shape, ReinitFunction reinit_f, OpSizeType size_type,
            std::vector<PatchOp<spacedim> *> input_ops = {}, uint n_dofs = 1)
    : dim_(dim), shape_(set_shape_vec(shape)), size_type_(size_type), input_ops_(input_ops),
      n_dofs_(n_dofs), reinit_func(reinit_f)
    {}

    /// Aligns shape_vec to 2 items (equal to matrix number of dimensions)
    std::vector<uint> set_shape_vec(std::initializer_list<uint> shape) const {
        std::vector<uint> shape_vec(shape);
        if (shape_vec.size() == 1) shape_vec.push_back(1);
        ASSERT_EQ(shape_vec.size(), 2);
        return shape_vec;
    }

    /**
     * Return number of operation components
     *
     * Value is computed from shape_ vector
     */
    inline uint n_comp() const {
        return shape_[0] * shape_[1];
    }

    /// Getter for dimension
    inline uint dim() const {
        return dim_;
    }

    /// Getter for size_type_
    OpSizeType size_type() const {
        return size_type_;
    }

    /// Getter for n_dofs_
    inline uint n_dofs() const {
        return n_dofs_;
    }

    /// Return pointer to operation of i_op index in input operation vector.
    inline PatchOp<spacedim> *input_ops(uint i_op) const {
        return input_ops_[i_op];
    }

    /// Getter for shape_
    inline const std::vector<uint> &shape() const {
        return shape_;
    }

    /**
     * Format shape to string
     *
     * Method is used in output development method.
     */
    inline std::string format_shape() const {
        std::stringstream ss;
        ss << shape_[0] << "x" << shape_[1];
        return ss.str();
    }

    /// Call reinit function on element table if function is defined
    inline void reinit_function(FMT_UNUSED std::vector<PatchOp<spacedim> *> &operations, IntTableArena &int_table) {
        reinit_func(this, int_table);
    }

    inline void allocate_result(size_t data_size, PatchArena &arena) {
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1]);
        for (uint i=0; i<n_comp(); ++i)
            result_(i) = ArenaVec<double>(data_size, arena);
    }

    inline void allocate_const_result(ArenaVec<double> &value_vec) {
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1]);
        for (uint i=0; i<n_comp(); ++i)
            result_(i) = value_vec;
    }

    /// Return map referenced result as Eigen::Vector
    Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> result_matrix() {
        return Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>>(result_.data(), shape_[0], shape_[1]);
    }

    /// Return map referenced result of DOF values as Eigen::Matrix
    Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> result_sub_matrix(uint i_dof) {
        ASSERT_LT(i_dof, n_dofs_);
        uint n_dof_comps = shape_[0] * shape_[1] / n_dofs_;
        return Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>>(result_.data()+i_dof*n_dof_comps, shape_[0], shape_[1] / n_dofs_);
    }

    /// Return map referenced result as Eigen::Vector
    Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> &raw_result() {
        return result_;
    }

    /// Same as previous but return const reference
    const Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> &raw_result() const {
        return result_;
    }


protected:
    uint dim_;                                    ///< Dimension
    std::vector<uint> shape_;                     ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> result_;    ///< Result matrix of operation
    OpSizeType size_type_;                         ///< Type of operation by size of vector (element, point or fixed size)
    std::vector<PatchOp<spacedim> *> input_ops_;  ///< Indices of operations in PatchPointValues::operations_ vector on which PatchOp is depended
    uint n_dofs_;                                 ///< Number of DOFs of FE operations (or 1 in case of element operations)

    ReinitFunction reinit_func;                   ///< Pointer to patch reinit function of element data table specialized by operation
};



/// Defines common functionality of reinit operations.
struct common_reinit {
	// empty base operation
	static inline void op_base(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // empty
    }

    template<unsigned int dim>
    static inline void elop_jac(PatchOp<3> * result_op) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto jac_value = result_op->result_matrix();
        auto coords_value = result_op->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<dim; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }

    template<unsigned int dim>
    static inline void elop_inv_jac(PatchOp<3> * result_op) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto inv_jac_value = result_op->result_matrix();
        auto jac_value = result_op->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<3, dim>(jac_value);
    }

    template<unsigned int dim>
    static inline void elop_jac_det(PatchOp<3> * result_op) {
        // result double, input matrix(spacedim, dim)
        auto jac_det_value = result_op->result_matrix();
        auto jac_value = result_op->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<3, dim>(jac_value).abs();
    }

    static inline void ptop_JxW(PatchOp<3> * result_op) {
        auto weights_value = result_op->input_ops(0)->result_matrix();
        auto jac_det_value = result_op->input_ops(1)->result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0,0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        auto jxw_value = result_op->result_matrix();
        jxw_value(0,0) = jxw_ovec.get_vec();
    }

    /// Common reinit function of vector symmetric gradient on bulk and side points
    static inline void ptop_vector_sym_grad(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        for (uint i_dof=0; i_dof<result_op->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = result_op->input_ops(0)->result_sub_matrix(i_dof);
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > sym_grad_dof = result_op->result_sub_matrix(i_dof);
            sym_grad_dof = 0.5 * (grad_vector_dof.transpose() + grad_vector_dof);
        }
    }
	/// Common reinit function of vector divergence on bulk and side points
    static inline void ptop_vector_divergence(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto divergence_value = result_op->result_matrix();

        for (uint i_dof=0; i_dof<result_op->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = result_op->input_ops(0)->result_sub_matrix(i_dof);
            divergence_value(i_dof) = grad_vector_dof(0,0) + grad_vector_dof(1,1) + grad_vector_dof(2,2);
        }
    }
};

/// Defines reinit operations on bulk points.
struct bulk_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_inv_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_inv_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_jac_det(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac_det<dim>(result_op);
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_JxW(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::ptop_JxW(result_op);
    }
    static inline void ptop_scalar_shape(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto ref_vec = result_op->input_ops(0)->result_matrix();
        auto result_vec = result_op->result_matrix();

        uint n_dofs = result_op->n_dofs();
        uint n_elem = result_vec(0).data_size() / ref_vec(0).data_size();

        ArenaVec<double> elem_vec(n_elem, result_op->raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> ref_ovec(n_dofs);
        for (uint i=0; i<n_dofs; ++i) {
            ref_ovec(i) = ArenaOVec<double>( ref_vec(i) );
        }

        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> result_ovec = elem_ovec * ref_ovec;
        for (uint i=0; i<n_dofs; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    static inline void ptop_vector_shape(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto ref_shape_vec = result_op->input_ops(0)->result_matrix();
        auto result_vec = result_op->result_matrix(); // result: shape 3x1

        uint n_dofs = result_op->n_dofs();
        uint n_elem = result_vec(0).data_size() / ref_shape_vec(0).data_size();

        ArenaVec<double> elem_vec(n_elem, result_op->raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_ovec(3, n_dofs);
        for (uint c=0; c<3*n_dofs; ++c) {
            ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(c));
        }

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = elem_ovec * ref_shape_ovec;
        for (uint c=0; c<3*n_dofs; ++c)
            result_vec(c) = result_ovec(c).get_vec();
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = op.input_ops(0).result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    static inline void ptop_vector_piola_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = op.input_ops(0).result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto inv_jac_vec = result_op->input_ops(0)->result_matrix();    // dim x spacedim=3
        auto ref_grads_vec = result_op->input_ops(1)->result_matrix();  // dim x n_dofs

        uint n_dofs = result_op->n_dofs();

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_grads_ovec(dim, n_dofs);
        for (uint i=0; i<dim*n_dofs; ++i) {
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = result_op->result_matrix();
        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<3*n_dofs; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto inv_jac_vec = result_op->input_ops(0)->result_matrix();   // dim x 3
        auto ref_grads_vec = result_op->input_ops(1)->result_matrix(); // dim x 3*n_dofs
        auto result_vec = result_op->result_matrix();                  // 3 x 3*n_dofs

        uint n_dofs = result_op->n_dofs();

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> ref_grads_ovec;
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i=0; i<3*dim; ++i) {
                ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i_dof*3*dim + i));
            }

            Eigen::Matrix<ArenaOVec<double>, 3, 3> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
            for (uint i=0; i<3; ++i) {
                for (uint j=0; j<3; ++j) {
                   result_vec(j,i+3*i_dof) = result_ovec(i,j).get_vec();
                }
            }
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
};



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_el_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_el_inv_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_inv_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        common_reinit::elop_jac<dim-1>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac_det(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac_det<dim-1>(result_op);
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_JxW(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::ptop_JxW(result_op);
    }
    template<unsigned int dim>
    static inline void ptop_normal_vec(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto normal_value = result_op->result_matrix();
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();
        normal_value = inv_jac_value.transpose() * RefElement<dim>::normal_vector_array( el_table(3) );

        ArenaVec<double> norm_vec( normal_value(0).data_size(), normal_value(0).arena() );
        Eigen::VectorXd A(3);
        for (uint i=0; i<normal_value(0).data_size(); ++i) {
            A(0) = normal_value(0)(i);
            A(1) = normal_value(1)(i);
            A(2) = normal_value(2)(i);
            norm_vec(i) = A.norm();
        }

        for (uint i=0; i<3; ++i) {
            normal_value(i) = normal_value(i) / norm_vec;
        }
    }
    static inline void ptop_scalar_shape(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto ref_vec = result_op->input_ops(0)->result_matrix();
        auto result_vec = result_op->result_matrix();

        uint n_dofs = result_op->n_dofs();
        uint n_sides = el_table(3).data_size();        // number of sides on patch
        uint n_patch_points = el_table(4).data_size(); // number of points on patch

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt) {
                result_vec(i_dof)(i_pt) = ref_vec(el_table(4)(i_pt), i_dof)(i_pt / n_sides);
            }
        }
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    static inline void ptop_vector_piola_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    static inline void ptop_vector_shape(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto ref_shape_vec = result_op->input_ops(0)->result_matrix();  // dim+1 x 3*n_dofs
        auto result_vec = result_op->result_matrix();                   // 3 x n_dofs

        uint n_dofs = result_op->n_dofs();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        for (uint c=0; c<3*n_dofs; c++)
        	result_vec(c) = ArenaVec<double>(n_patch_points, result_vec(c).arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                for (uint c=0; c<3; c++)
                    result_vec(c,i_dof)(i_pt) = ref_shape_vec(el_table(4)(i_pt),3*i_dof+c)(i_pt / n_sides);
        }
    }

template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto ref_shape_grads = result_op->input_ops(1)->result_matrix();  // dim+1 x dim*n_dofs
        auto grad_scalar_shape_value = result_op->result_matrix();  // Result vector: 3 x n_dofs

        uint n_dofs = result_op->n_dofs();
        uint n_points = ref_shape_grads(0).data_size();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_patch_points, inv_jac_value(i).arena() );
        	for (uint j=0; j<n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_grads_expd(dim, n_dofs);
        for (uint i=0; i<dim*n_dofs; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, inv_jac_value(0).arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = i_pt * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<dim; ++i_c) {
                        ref_shape_grads_expd(i_c, i_dof)(i_begin + i_sd) = ref_shape_grads(el_table(3)(i_sd), i_dof*dim+i_c)(i_pt);
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto result_vec = result_op->raw_result();                        // 3 x 3*n_dofs
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();    // dim x 3
        auto ref_vector_grad = result_op->input_ops(1)->result_matrix();  // n_sides*dim x 3*n_dofs

        uint n_dofs = result_op->n_dofs();
        uint n_points = ref_vector_grad(0).data_size();
        uint n_patch_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_patch_points, inv_jac_value(i).arena() );
        	for (uint j=0; j<n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_patch_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, dim, 3> ref_shape_grads_expd;
        for (uint i=0; i<3*dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, inv_jac_value(0).arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {

            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = i_pt * n_patch_sides;
                for (uint i_sd=0; i_sd<n_patch_sides; ++i_sd) {
                    for (uint i_dim=0; i_dim<dim; ++i_dim) {
                        for (uint i_c=0; i_c<3; ++i_c) {
                            ref_shape_grads_expd(i_dim, i_c)(i_begin + i_sd) = ref_vector_grad(el_table(3)(i_sd)*dim+i_dim, 3*i_dof+i_c)(i_pt);
                        }
                    }
                }
            }

            // computes operation result
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > res_submat = result_op->result_sub_matrix(i_dof);
            res_submat = (inv_jac_expd_value.transpose() * ref_shape_grads_expd).transpose();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
};


// template specialization
template<>
inline void side_reinit::elop_sd_jac<1>(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
}

template<>
inline void side_reinit::elop_sd_jac_det<1>(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
    auto result_vec = result_op->result_matrix();
    for (uint i=0;i<result_vec(0).data_size(); ++i) {
        result_vec(0,0)(i) = 1.0;
    }
}


#endif /* OP_FUNCTION_HH_ */
