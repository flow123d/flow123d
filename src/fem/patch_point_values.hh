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


enum PointType {
    bulk_point,
    side_point
};


template<unsigned int spacedim = 3>
class PatchPointValues
{
public:
    /// Default constructor, set invalid dim (uninitialized object)
    PatchPointValues()
    : dim_(10) {}

    /// Initialize object, set number of columns (quantities) in tables
    void initialize(uint dim, PointType point_type, uint point_cols, uint int_cols) {
        ASSERT_EQ(dim_, 10).error("Multiple initialization!\n");
        ASSERT(dim <= 3)(dim).error("Dimension must be 0, 1, 2 or 3!\n");

        this->reset();
        dim_ = dim;
    	point_type_ = point_type;

    	point_vals_.resize(point_cols);
    	int_vals_.resize(int_cols);
    	el_vals_.resize( (dim_+1) * spacedim );
    }

    /// Reset number of rows (points)
    void reset() {
    	ASSERT(dim_ <= 3).error("Uninitialized PatchPointValues object!\n");
        n_points_ = 0;
        n_elems_ = 0;
    }

private:
    TableDbl point_vals_;
    TableInt int_vals_;
    TableDbl el_vals_;

    uint dim_;
    PointType point_type_;

    uint n_points_;
    uint n_elems_;
};


/**
 * Base class of all FE operations.
 */
class FeOp {
public:
    FeOp(bool enabled, bool input_column, PatchPointValues &point_vals)
    : enabled_(enabled), input_column_(input_column), point_vals_(point_vals)
    {}

    inline Scalar scalar_val(uint point_idx) const {
        return point_vals_[input_column_][point_vals];
    }

    inline Vector vector_val(uint point_idx) const {
        Vector val;
        for (uint i=0; i<3; ++i)
            val(i) = point_vals_[input_column_+i][point_vals];
        return val;
    }

    inline Tensor tensor_val(uint point_idx) const {
        Tensor val;
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = point_vals_[input_column_+3*i+j][point_vals];
        return val;
    }


private:
    bool enabled_;                  ///< if the operation is active, column will be updated
    bool input_column_;             ///< first column to scalar, vector or matrix inputs
    PatchPointValues &point_vals_;  ///< Reference to data table
};

#endif /* PATCH_POINT_VALUES_HH_ */
