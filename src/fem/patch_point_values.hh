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




/**
 * @brief Class for storing FE data of quadrature points.
 *
 * Store data of bulk or side quadrature points of one dimension.
 */
template<unsigned int spacedim = 3>
class PatchPointValues
{
public:
    /**
     * Constructor
     *
     * Set dimension
     */
    PatchPointValues(uint dim);

    /**
     * Initialize object, set number of columns (quantities) in tables.
     *
     * Number of columns of int_vals_ table is passed by argument, number of columns
     * of other tables is given by n_columns_ value.
     */
    void initialize(uint int_cols);

    /// Reset number of rows (points and elements)
    inline void reset() {
        n_points_ = 0;
        n_elems_ = 0;
    }

    /// Getter of n_columns_
    inline uint n_columns() const {
        return n_columns_;
    }

    /// Adds the number of columns equal to n_added, returns index of first of them
    /// Temporary method allow acces to old structure PatchFEValues::DimPatchFEValues
    inline uint add_columns(uint n_added) {
        uint old_size = n_columns_;
        n_columns_ += n_added;
        return old_size;
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(arma::mat coords, uint element_patch_idx);

    /// Add accessor to operations_ vector
    ElOp<spacedim> &add_accessor(ElOp<spacedim> op_accessor);

protected:
    /**
     * Store data of bulk or side quadrature points of one dimension
     *
     * Number of columns is given by n_columns_, number of used rows by n_points_.
     */
    TableDbl point_vals_;
    /**
     * Hold integer values of quadrature points of previous table.
     *
     * Table contains 2 columns for bulk point table (index of point on patch, element_idx)
     * and 3 columns for side point table (index of point on patch, element_idx, side_idx).
     * Number of used rows is given by n_points_.
     */
    TableInt int_vals_;
    /**
     * Store data of elements. Elements are given by quadrature points stored in point_vals_
     *
     * Number of columns is given by n_columns_, number of used rows by n_elems_.
     */
    TableDbl el_vals_;

    /// Vector of all defined operations
    std::vector<ElOp<spacedim>> operations_;

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
 * @brief Class represents FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
    /// Type for conciseness
    typedef void (*ReinitFunction)(std::vector<ElOp<spacedim> *> &, TableDbl &);

    /// Constructor
    ElOp(uint dim, std::initializer_list<uint> shape, uint result_col, ReinitFunction r_func, ElOp<spacedim> *input_op = nullptr)
    : dim_(dim), shape_(shape), result_col_(result_col), reinit_func(r_func)
    {
        if (input_op != nullptr) {
            uint first_col = input_op->result_col();
            for (uint i=0; i<input_op->n_comp(); ++i)
                input_column_.push_back(i+first_col);
        }
    }

    /// Number of components computed from shape_ vector
    inline uint n_comp() const {
        if (shape_.size() == 1) return shape_[0];
        else return shape_[0] * shape_[1];
    }

    /// Getter of result_col_
    inline uint result_col() const {
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
    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    uint result_col_;                         ///< First column to scalar, vector or matrix result
    std::vector<uint> input_column_;          ///< Vector of column on which ElOp is depended

    /// Pointer to patch reinit function specialized by operation
    ReinitFunction reinit_func;
};


namespace FeBulk {

/// enum of bulk operation
enum BulkOps
{
	opCoords,
	opElCoords,
	opJac,
	opJacDet
};

/// Bulk data specialization, order of item in operations_ vector corresponds to the BulkOps enum
template<unsigned int spacedim = 3>
class PatchPointValues : public ::PatchPointValues<spacedim> {
public:
    /// Constructor fill operations_ vector.
    PatchPointValues(uint dim);
};

} // closing namespace FeBulk

/// Defines reinit operations on bulk points.
struct bulk_ops {
    static inline void reinit_elop_coords(std::vector<ElOp<3> *> &operations, TableDbl &op_results) {
        std::cout << operations.size() << " - " << op_results[0][0] << std::endl;
    }
    static inline void reinit_ptop_coords(std::vector<ElOp<3> *> &operations, TableDbl &op_results) {
        std::cout << operations.size() << " - " << op_results[0][0] << std::endl;
    }
    static inline void reinit_elop_jac(std::vector<ElOp<3> *> &operations, TableDbl &op_results) {
        std::cout << operations.size() << " - " << op_results[0][0] << std::endl;
    }
    static inline void reinit_elop_jac_det(std::vector<ElOp<3> *> &operations, TableDbl &op_results) {
        std::cout << operations.size() << " - " << op_results[0][0] << std::endl;
    }
};



namespace FeSide {

template<unsigned int spacedim = 3>
class PatchPointValues : public ::PatchPointValues<spacedim> {
public:
    /// Default constructor
    PatchPointValues(uint dim);
};

} // closing namespace FeSide


#endif /* PATCH_POINT_VALUES_HH_ */
