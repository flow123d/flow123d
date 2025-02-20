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

#include "fields/eval_subset.hh"

template<unsigned int spacedim> class PatchPointValues;
template<unsigned int spacedim> class PatchOp;

using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;



template <class ValueType>
class ElQ {
public:
    /// Forbidden default constructor
    ElQ() = delete;

    /// Constructor
    ElQ(PatchPointValues<3> *patch_point_vals, unsigned int op_idx)
    : patch_point_vals_(patch_point_vals), op_idx_(op_idx), patch_op_(nullptr) {}

    /// Constructor
    ElQ(PatchOp<3> *op)
    : patch_point_vals_(nullptr), op_idx_(0), patch_op_(op) {}

    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

private:
    PatchPointValues<3> *patch_point_vals_; ///< Pointer to PatchPointValues - temporary
    unsigned int op_idx_;                   ///< Index of operation in patch_point_vals_.operations vector - temporary
    PatchOp<3> *patch_op_;                  ///< Pointer to operation
};


template <class ValueType>
class FeQ {
public:
    /// Forbidden default constructor
    FeQ() = delete;

    // Class similar to current FeView - obsolete
    FeQ(PatchPointValues<3> *patch_point_vals, bool is_bulk, unsigned int op_idx)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(op_idx), patch_op_bulk_(nullptr), patch_op_side_(nullptr), i_shape_fn_idx_(0) {
        if (is_bulk) patch_point_vals_bulk_ = patch_point_vals;
        else patch_point_vals_side_ = patch_point_vals;
    }

    /// Constructor used only in FeQArray::shape() - obsolete
    FeQ(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side,
            unsigned int op_idx, unsigned int i_shape_fn_idx)
    : patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side),
      patch_op_bulk_(nullptr), patch_op_side_(nullptr), op_idx_(op_idx), i_shape_fn_idx_(i_shape_fn_idx) {}

    // Class similar to current FeView
    FeQ(PatchOp<3> *patch_op, bool is_bulk)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(0), patch_op_bulk_(nullptr), patch_op_side_(nullptr), i_shape_fn_idx_(0) {
        if (is_bulk) patch_op_bulk_ = patch_op;
        else patch_op_side_ = patch_op;
    }

    /// Constructor used only in FeQArray::shape()
    FeQ(PatchOp<3> *patch_op_bulk, PatchOp<3> *patch_op_side, unsigned int i_shape_fn_idx)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(0),
      patch_op_bulk_(patch_op_bulk), patch_op_side_(patch_op_side), i_shape_fn_idx_(i_shape_fn_idx) {}


    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

private:
    PatchPointValues<3> *patch_point_vals_bulk_; ///< Pointer to bulk PatchPointValues - obsolete
    PatchPointValues<3> *patch_point_vals_side_; ///< Pointer to side PatchPointValues - obsolete
    unsigned int op_idx_;                        ///< Index of operation in patch_point_vals_.operations vector - obsolete
    PatchOp<3> *patch_op_bulk_;                  ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                  ///< Pointer to side operation
    unsigned int i_shape_fn_idx_;                ///< Index of shape function
};


template <class ValueType>
class FeQArray {
public:
    /// Forbidden default constructor
    FeQArray() = delete;

    // Class similar to current FeView
    FeQArray(PatchPointValues<3> *patch_point_vals, bool is_bulk, unsigned int op_idx, unsigned int n_dofs)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(op_idx), patch_op_bulk_(nullptr), patch_op_side_(nullptr), n_dofs_(n_dofs) {
        ASSERT_GT(n_dofs, 0).error("Invalid number of DOFs.\n");

        if (is_bulk) patch_point_vals_bulk_ = patch_point_vals;
        else patch_point_vals_side_ = patch_point_vals;
    }

    FeQArray(PatchOp<3> *patch_op, bool is_bulk)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(0), patch_op_bulk_(nullptr),
      patch_op_side_(nullptr), n_dofs_(patch_op->n_dofs()) {
        ASSERT_GT(n_dofs_, 0).error("Invalid number of DOFs.\n");

        if (is_bulk) patch_op_bulk_ = patch_op;
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
    PatchPointValues<3> *patch_point_vals_bulk_; ///< Reference to bulk PatchPointValues - obsolete
    PatchPointValues<3> *patch_point_vals_side_; ///< Reference to side PatchPointValues - obsolete
    unsigned int op_idx_;                        ///< Index of operation in patch_point_vals_.operations vector - obsolete
    PatchOp<3> *patch_op_bulk_;                  ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                  ///< Pointer to side operation
    unsigned int n_dofs_;                        ///< Number of DOFs
};


template <class ValueType>
class FeQJoin {
public:
    /// Default constructor
    FeQJoin()
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr),
      patch_op_bulk_(nullptr), patch_op_side_(nullptr), patch_op_zero_bulk_(nullptr), patch_op_zero_side_(nullptr)
    {}

    /**
     * Constructor
     *
     * @param patch_point_vals_bulk  Pointer to PatchPointValues bulk object.
     * @param patch_point_vals_side  Pointer to PatchPointValues side object.
     * @param begin                  Index of the first component of the bulk Quantity.
     * @param begin_side             Index of the first component of the side Quantity.
     * @param n_dofs_bulk            Number of DOFs of bulk (lower-dim) element.
     * @param n_dofs_side            Number of DOFs of side (higher-dim) element.
     */
    FeQJoin(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side, unsigned int n_dofs_bulk,
            unsigned int n_dofs_side, unsigned int op_idx_bulk, unsigned int op_idx_side)
    : patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side),
      patch_op_bulk_(nullptr), patch_op_side_(nullptr), patch_op_zero_bulk_(nullptr), patch_op_zero_side_(nullptr),
      n_dofs_high_(n_dofs_side), n_dofs_low_(n_dofs_bulk), op_idx_bulk_(op_idx_bulk), op_idx_side_(op_idx_side) {}

    /**
     * Constructor
     *
     * @param patch_op_bulk  Pointer to PatchOP bulk object.
     * @param patch_op_side  Pointer to PatchOp side object.
     */
    FeQJoin(PatchOp<3> *patch_op_bulk, PatchOp<3> *patch_op_side, PatchOp<3> *patch_op_zero_bulk, PatchOp<3> *patch_op_zero_side)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr),
      patch_op_bulk_(patch_op_bulk), patch_op_side_(patch_op_side), patch_op_zero_bulk_(patch_op_zero_bulk), patch_op_zero_side_(patch_op_zero_side),
      n_dofs_high_(patch_op_side->n_dofs()), n_dofs_low_(patch_op_bulk->n_dofs()), op_idx_bulk_(0), op_idx_side_(0) {}


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
    PatchPointValues<3> *patch_point_vals_bulk_;  ///< Pointer to bulk PatchPointValues - obsolete
    PatchPointValues<3> *patch_point_vals_side_;  ///< Pointer to side PatchPointValues - obsolete
    PatchOp<3> *patch_op_bulk_;                   ///< Pointer to bulk operation
    PatchOp<3> *patch_op_side_;                   ///< Pointer to side operation
    PatchOp<3> *patch_op_zero_bulk_;              ///< Pointer to bulk zero operation
    PatchOp<3> *patch_op_zero_side_;              ///< Pointer to side zero operation
    unsigned int n_dofs_high_;                    ///< Number of DOFs on high-dim element
    unsigned int n_dofs_low_;                     ///< Number of DOFs on low-dim element
    unsigned int op_idx_bulk_;                    ///< Index of operation in patch_point_vals_bulk_.operations vector - obsolete
    unsigned int op_idx_side_;                    ///< Index of operation in patch_point_vals_side_.operations vector - obsolete
};


#endif /* OP_ACCESSORS_HH_ */
