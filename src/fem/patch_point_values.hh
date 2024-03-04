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
 * @file    patch_point_values.hh
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#ifndef PATCH_POINT_VALUES_HH_
#define PATCH_POINT_VALUES_HH_

#include "fem/eigen_tools.hh"
#include "fem/dh_cell_accessor.hh"


template<unsigned int spacedim> class PatchFEValues;
template<unsigned int spacedim> class ElOp;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;


static const uint INVALID_COLUMN = -1;



template<unsigned int spacedim = 3>
class PatchPointValues
{
public:
    /// Default constructor
    PatchPointValues(uint dim);

    /// Initialize object, set number of columns (quantities) in tables
    void initialize(uint int_cols);

    /// Reset number of rows (points)
    inline void reset() {
        n_points_ = 0;
        n_elems_ = 0;
    }

    /// Getter of n_columns_
    inline uint n_columns() const {
        return n_columns_;
    }

    /// Adds the number of columns equal to n_added, returns index of first of them
    inline uint add_columns(uint n_added) {
        uint old_size = n_columns_;
        n_columns_ += n_added;
        return old_size;
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx);

    ElOp<spacedim> *add_accessor(ElOp<spacedim> *op_accessor);

protected:
    TableDbl point_vals_;
    TableInt int_vals_;
    TableDbl el_vals_;

    std::vector<ElOp<spacedim> *> operation_columns_;

    uint dim_;                        ///< Dimension
    uint n_columns_;                  ///< Number of columns of \p point_vals table
    uint n_points_;                   ///< Number of points in patch
    uint n_elems_;                    ///< Number of elements in patch

    std::vector<uint> elements_map_;  ///< Map of element patch indices to el_vals_ table
    std::vector<uint> points_map_;    ///< Map of point patch indices  to point_vals_ and int_vals_ tables

    friend class PatchFEValues<spacedim>;
    friend class ElOp<spacedim>;
};


/**
 * Base class of all FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
    /// Constructor
    ElOp(uint dim, std::initializer_list<uint> shape)
    : dim_(dim), shape_(shape), result_col_(INVALID_COLUMN), input_column_(INVALID_COLUMN)
    {}

    /// Reinit data on all patch points.
    virtual void reinit_data() =0;

    /// Number of components computed from shape_ vector
    inline uint n_comp() const {
        if (shape_.size() == 1) return shape_[0];
        else return shape_[0] * shape_[1];
    }

    /**
     * Register operation columns if result_col_ is not set.
     *
     * Return index of first column it point values table.
     */
    inline uint register_columns(PatchPointValues<spacedim> &point_vals) {
        ASSERT_EQ(dim_, point_vals.dim_);
        if (result_col_ == INVALID_COLUMN) {
            this->check_op_dependency(point_vals);
            result_col_ = point_vals.add_columns(this->n_comp());
        }
        return result_col_;
    }

//    inline Scalar scalar_val(uint point_idx) const {
//        return point_vals_->point_vals_[input_column_][point_idx];
//    }
//
//    inline Vector vector_val(uint point_idx) const {
//        Vector val;
//        for (uint i=0; i<3; ++i)
//            val(i) = point_vals_->point_vals_[input_column_+i][point_idx];
//        return val;
//    }
//
//    inline Tensor tensor_val(uint point_idx) const {
//        Tensor val;
//        for (uint i=0; i<3; ++i)
//            for (uint j=0; j<3; ++j)
//                val(i,j) = point_vals_->point_vals_[input_column_+3*i+j][point_idx];
//        return val;
//    }


protected:
    /// Check and register dependency of operations, can be implemented in descendants.
    virtual void check_op_dependency(FMT_UNUSED PatchPointValues<spacedim> &point_vals) {}

    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    uint result_col_;                         ///< Result column.
    uint input_column_;                       ///< First column to scalar, vector or matrix inputs
};


namespace FeBulk {

enum BulkOps
{
	opCoords,
	opElCoords,
	opJac,
	opJacDet
};

template<unsigned int spacedim = 3>
class PatchPointValues : public ::PatchPointValues<spacedim> {
public:
    /// Default constructor
    PatchPointValues(uint dim);
};


template<unsigned int spacedim = 3>
class OpCoords : public ElOp<spacedim> {
public:
	/// Constructor
	OpCoords(uint dim)
    : ElOp<spacedim>(dim, {spacedim})
    {}

	/// Implement ElOp::reinit_data
	void reinit_data() override;
};


template<unsigned int spacedim = 3>
class OpElCoords : public ElOp<spacedim> {
public:
    /// Constructor
    OpElCoords(uint dim)
    : ElOp<spacedim>(dim, {spacedim, dim+1})
    {}

	/// Implement ElOp::reinit_data
	void reinit_data() override;
};


template<unsigned int spacedim = 3>
class OpJac : public ElOp<spacedim> {
public:
    /// Constructor
    OpJac(uint dim, ElOp<spacedim> *coords_op)
    : ElOp<spacedim>(dim, {spacedim, dim}), coords_operator_(coords_op)
    {}

	/// Implement ElOp::reinit_data
	void reinit_data() override;

	void check_op_dependency(::PatchPointValues<spacedim> &point_vals) override;
private:
	ElOp<spacedim> *coords_operator_;
};


template<unsigned int spacedim = 3>
class OpJacDet : public ElOp<spacedim> {
public:
    /// Constructor
    OpJacDet(uint dim, ElOp<spacedim> *jac_op)
    : ElOp<spacedim>(dim, {1}), jac_operator_(jac_op)
    {}

    /// Implement ElOp::reinit_data
    void reinit_data() override;

	void check_op_dependency(::PatchPointValues<spacedim> &point_vals) override;

private:
	ElOp<spacedim> *jac_operator_;
};

} // closing namespace FeBulk


namespace FeSide {

template<unsigned int spacedim = 3>
class PatchPointValues : public ::PatchPointValues<spacedim> {
public:
    /// Default constructor
    PatchPointValues(uint dim);
};

} // closing namespace FeSide


#endif /* PATCH_POINT_VALUES_HH_ */
