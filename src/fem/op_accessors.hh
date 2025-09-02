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
 * @file    op_accessors.hh
 * @brief   Declares accessors to FE operations.
 * @author  David Flanderka
 */

#ifndef OP_ACCESSORS_HH_
#define OP_ACCESSORS_HH_

#include "fem/integral_points.hh"
#include "fem/patch_op.hh"

using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;



template <class ValueType>
class ElQ {
public:
    /**
     * Default constructor.
     *
     * Used only in template specialization of CoplingIntegral. DO NOT USE in other cases.
     */
    ElQ()
    : patch_op_(nullptr) {}

    /// Constructor
    ElQ(PatchOp<3> *op)
    : patch_op_(op) {}

    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

private:
    PatchOp<3> *patch_op_;                  ///< Pointer to operation
};


template <class ValueType>
class FeQ {
public:
    /**
     * Default constructor.
     *
     * Used only in template specialization of CoplingIntegral. DO NOT USE in other cases.
     */
    FeQ()
    : patch_op_bulk_(nullptr), patch_op_side_(nullptr), i_shape_fn_idx_(0) {}

    // Class similar to current FeView
    FeQ(PatchOp<3> *patch_op)
    : patch_op_bulk_(nullptr), patch_op_side_(nullptr), i_shape_fn_idx_(0) {
        if (patch_op->domain()==0) patch_op_bulk_ = patch_op;
        else patch_op_side_ = patch_op;
    }

    /// Constructor used only in FeQArray::shape()
    FeQ(PatchOp<3> *patch_op_bulk, PatchOp<3> *patch_op_side, unsigned int i_shape_fn_idx)
    : patch_op_bulk_(patch_op_bulk), patch_op_side_(patch_op_side), i_shape_fn_idx_(i_shape_fn_idx) {}


    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

private:
    PatchOp<3> *patch_op_bulk_;                  ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                  ///< Pointer to side operation
    unsigned int i_shape_fn_idx_;                ///< Index of shape function
};


template <class ValueType>
class FeQArray {
public:
    /// Default constructor
    FeQArray()
    : patch_op_bulk_(nullptr), patch_op_side_(nullptr), n_dofs_(1) {}

    // Class similar to current FeView
    FeQArray(PatchOp<3> *patch_op)
    : patch_op_bulk_(nullptr), patch_op_side_(nullptr), n_dofs_(patch_op->n_dofs()) {
        ASSERT_GT(n_dofs_, 0).error("Invalid number of DOFs.\n");

        if (patch_op->domain()==0) patch_op_bulk_ = patch_op;
        else patch_op_side_ = patch_op;
    }


    FeQ<ValueType> shape(unsigned int i_shape_fn_idx) const {
        ASSERT_LT(i_shape_fn_idx, n_dofs_);
        return FeQ<ValueType>(patch_op_bulk_, patch_op_side_, i_shape_fn_idx);
    }

    /// Return number of DOFs
    inline unsigned int n_dofs() const {
        return n_dofs_;
    }

private:
    PatchOp<3> *patch_op_bulk_;                  ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                  ///< Pointer to side operation
    unsigned int n_dofs_;                        ///< Number of DOFs
};


template <class ValueType>
class FeQJoin {
public:
    /// Default constructor
    FeQJoin()
    : patch_op_bulk_(nullptr), patch_op_side_(nullptr), patch_op_zero_bulk_(nullptr), patch_op_zero_side_(nullptr)
    {}

    /**
     * Constructor
     *
     * @param patch_op_bulk       Pointer to PatchOP bulk object.
     * @param patch_op_side       Pointer to PatchOp side object.
     * @param patch_op_zero_bulk  Pointer to zero PatchOP bulk object.
     * @param patch_op_zero_side  Pointer to zero PatchOp side object.
     */
    FeQJoin(PatchOp<3> *patch_op_bulk, PatchOp<3> *patch_op_side, PatchOp<3> *patch_op_zero_bulk, PatchOp<3> *patch_op_zero_side)
    : patch_op_bulk_(patch_op_bulk), patch_op_side_(patch_op_side), patch_op_zero_bulk_(patch_op_zero_bulk), patch_op_zero_side_(patch_op_zero_side),
      n_dofs_high_(patch_op_side->n_dofs()), n_dofs_low_(patch_op_bulk->n_dofs()) {}


    inline unsigned int n_dofs_low() const {
        return n_dofs_low_;
    }

    inline unsigned int n_dofs_high() const {
        return n_dofs_high_;
    }

    inline unsigned int n_dofs_both() const {
        return n_dofs_high_ + n_dofs_low_;
    }

//    /// Return local index of DOF (on low / high-dim) - should be private method
//    inline unsigned int local_idx(unsigned int i_join_idx) const {
//        if (this->is_high_dim(i_join_idx)) return (i_join_idx - n_dofs_low());
//        else return i_join_idx;
//    }

    inline bool is_high_dim(unsigned int i_join_idx) const {
        return (i_join_idx >= n_dofs_low());
    }

    FeQ<ValueType> shape(unsigned int i_join_idx) const {
        ASSERT_LT(i_join_idx, n_dofs_both());

        /*
         * Set zero bulk PatchFeValues for side DOFs and zero side PatchFeValues for bulk DOFs
         *
         * TODO After implementation of dependencies:
         *      1) Implement FeQ::vec() getter (experimental method returned entire data vector)
         *      2) Test difference of vectors
         */
        if (this->is_high_dim(i_join_idx))
            return FeQ<ValueType>(patch_op_zero_side_, patch_op_side_, i_join_idx - n_dofs_low());
        else
            return FeQ<ValueType>(patch_op_bulk_, patch_op_zero_bulk_, i_join_idx);
    }


private:
    // attributes:
    PatchOp<3> *patch_op_bulk_;                   ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                   ///< Pointer to side operation
    PatchOp<3> *patch_op_zero_bulk_;              ///< Pointer to bulk zero operation
    PatchOp<3> *patch_op_zero_side_;              ///< Pointer to side zero operation
    unsigned int n_dofs_high_;                    ///< Number of DOFs on high-dim element
    unsigned int n_dofs_low_;                     ///< Number of DOFs on low-dim element
};


#endif /* OP_ACCESSORS_HH_ */
