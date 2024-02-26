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


template<unsigned int spacedim> class ElOp;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;



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

    ElOp<spacedim> &add_accessor(ElOp<spacedim> op_accessor);

protected:
    TableDbl point_vals_;
    TableInt int_vals_;
    TableDbl el_vals_;

    std::vector<ElOp<spacedim>> operation_columns_;

    uint dim_;

    uint n_columns_;       ///< Number of columns of \p point_vals table
    uint n_points_;
    uint n_elems_;

    friend class ElOp<spacedim>;
};


/**
 * Base class of all FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
	/// Constructor
	ElOp(uint dim, std::initializer_list<uint> shape, PatchPointValues<spacedim> &point_vals)
    : dim_(dim), shape_(shape), result_col_(-1), input_column_(-1), point_vals_(point_vals)
    {}

	/// Reinit data on all patch points.
	virtual void reinit_data() {
		ASSERT_PERMANENT(false).error("Method must be implemented in descendants");
	}

    inline Scalar scalar_val(uint point_idx) const {
        return point_vals_.point_vals_[input_column_][point_idx];
    }

    inline Vector vector_val(uint point_idx) const {
        Vector val;
        for (uint i=0; i<3; ++i)
            val(i) = point_vals_.point_vals_[input_column_+i][point_idx];
        return val;
    }

    inline Tensor tensor_val(uint point_idx) const {
        Tensor val;
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = point_vals_.point_vals_[input_column_+3*i+j][point_idx];
        return val;
    }


protected:
    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    uint result_col_;                         ///< Result column.
    uint input_column_;                       ///< First column to scalar, vector or matrix inputs
    PatchPointValues<spacedim> &point_vals_;  ///< Reference to data table
};


namespace FeBulk {

template<unsigned int spacedim = 3>
class PatchPointValues : public ::PatchPointValues<spacedim> {
public:
    /// Default constructor
    PatchPointValues(uint dim);
};


template<unsigned int spacedim = 3>
class OpJac : public ElOp<spacedim> {
public:
	/// Constructor
	OpJac(uint dim, PatchPointValues<spacedim> &point_vals)
    : ElOp<spacedim>(dim, {spacedim, dim}, point_vals)
    {}


	/// Implement ElOp::reinit_data
	void reinit_data() override;
};


template<unsigned int spacedim = 3>
class OpJacDet : public ElOp<spacedim> {
public:
	/// Constructor
	OpJacDet(uint dim, PatchPointValues<spacedim> &point_vals, ElOp<spacedim> &jac_op)
    : ElOp<spacedim>(dim, {1}, point_vals), jac_operator_(jac_op)
    {}

    /// Implement ElOp::reinit_data
    void reinit_data() override;

private:
	ElOp<spacedim> &jac_operator_;
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
