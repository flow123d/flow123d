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

#include <Eigen/Dense>

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

    /// Getter of n_elems_
    inline uint n_elems() const {
        return n_elems_;
    }

    /// Getter of n_points_
    inline uint n_points() const {
        return n_points_;
    }

    /// Adds the number of columns equal to n_added, returns index of first of them
    /// Temporary method allow acces to old structure PatchFEValues::DimPatchFEValues
    inline uint add_columns(uint n_added) {
        uint old_size = n_columns_;
        n_columns_ += n_added;
        return old_size;
    }

    /// Register element to el_vals_ table
    uint register_element(arma::mat coords, uint element_patch_idx);

    /// Register side to el_vals_ table
    uint register_side(arma::mat coords);

    /// Register bulk point to point_vals_ and int_vals_ table
    uint register_bulk_point(uint elem_table_row, uint value_patch_idx, uint elem_idx);

    /// Register side point to point_vals_ and int_vals_ table
    uint register_side_point(uint elem_table_row, uint value_patch_idx, uint elem_idx, uint side_idx);

    /// Add accessor to operations_ vector
    ElOp<spacedim> &add_accessor(ElOp<spacedim> op_accessor);

    /// Reinit data.
    void reinit_patch() {
        // precompute data on el_vals_ table
        for (uint i=0; i<operations_.size(); ++i)
            operations_[i].reinit_function(operations_, el_vals_);

        // copy data from el_vals_ to point_vals_
        std::vector<uint> copied; // list of columns that will be copied
        for (uint i_op=0; i_op<operations_.size(); ++i_op)
            if (operations_[i_op].copy_vals()) {
            	uint res_col = operations_[i_op].result_col();
            	for (uint i_col=res_col; i_col<res_col+operations_[i_op].n_comp(); ++i_col)
            		copied.push_back(i_col);
            }
        for (uint i_pt=0; i_pt<n_points_; ++i_pt) {
            uint el_table_idx = int_vals_(1)(i_pt);
            for (uint i_q=0; i_q<copied.size(); ++i_q)
                point_vals_(copied[i_q])(i_pt) = el_vals_(copied[i_q])(el_table_idx);
        }

        // call reinit on point vals
    }

    /// Temporary development method
    void print(bool points, bool ints, bool elems) const {
        std::cout << "** Dimension: " << dim_ << std::endl;
        if (points) {
            std::cout << "Point vals: " << point_vals_.rows() << " - " << point_vals_.cols() << std::endl;
	        for (uint i_row=0; i_row<n_points_; ++i_row) {
                for (uint i_col=0; i_col<n_columns_; ++i_col)
                	std::cout << point_vals_(i_col)(i_row) << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (ints) {
            std::cout << "Int vals: " << int_vals_.rows() << " - " << int_vals_.cols() << std::endl;
	        for (uint i_row=0; i_row<n_points_; ++i_row) {
                for (uint i_col=0; i_col<3; ++i_col)
                	std::cout << int_vals_(i_col)(i_row) << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (elems) {
            std::cout << "El vals: " << el_vals_.rows() << " - " << el_vals_.cols() << std::endl;
	        for (uint i_row=0; i_row<n_elems_; ++i_row) {
                for (uint i_col=0; i_col<n_columns_; ++i_col)
                	std::cout << el_vals_(i_col)(i_row) << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "*****************" << std::endl;
    }

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
     * Table contains following columns:
     *  0: Index of quadrature point on patch
     *  1: Row of element/side in el_vals_ table
     *  2: Element idx in Mesh
     *  3: Index of side in element (column is allocated only for side point table)
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
    typedef void (*ReinitFunction)(std::vector<ElOp<spacedim>> &, TableDbl &);

    /// Constructor
    ElOp(uint dim, std::initializer_list<uint> shape, uint result_col, bool copy_vals, ReinitFunction r_func = nullptr, ElOp<spacedim> *input_op = nullptr)
    : dim_(dim), shape_(shape), result_col_(result_col), copy_vals_(copy_vals), reinit_func(r_func)
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

    /// Getter of dimension
    inline uint dim() const {
        return dim_;
    }

    /// Getter of result_col_
    inline uint result_col() const {
        return result_col_;
    }

    /// Getter of copy_vals_
    inline bool copy_vals() const {
        return copy_vals_;
    }

    inline void reinit_function(std::vector<ElOp<spacedim>> &operations, TableDbl &data_table) {
    	if (reinit_func != nullptr) reinit_func(operations, data_table);
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
    bool copy_vals_;                          ///< Flag marks if values of result columns are copied from el_vals to point_vals table

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
struct bulk_reinit {
	// element operations
	static inline void elop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void elop_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        uint dim = operations[FeBulk::BulkOps::opJac].dim();
        uint result_begin_col = operations[FeBulk::BulkOps::opJac].result_col();
        uint input_begin_col = operations[FeBulk::BulkOps::opElCoords].result_col();
        switch (dim) {
            case 1: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 1>> result_mat(op_results.data() + result_begin_col, 3, 1);
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> input_mat(op_results.data() + input_begin_col, 3, 2);
                result_mat = eigen_tools::jacobian<3,1>(input_mat);
                break;
            }
            case 2: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> result_mat(op_results.data() + result_begin_col, 3, 2);
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 3>> input_mat(op_results.data() + input_begin_col, 3, 3);
                result_mat = eigen_tools::jacobian<3,2>(input_mat);
                break;
            }
            case 3: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 3>> result_mat(op_results.data() + result_begin_col, 3, 3);
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 4>> input_mat(op_results.data() + input_begin_col, 3, 4);
                result_mat = eigen_tools::jacobian<3,3>(input_mat);
                break;
            }
        }
    }
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results) {
        // result double, input matrix(spacedim, dim)
        uint dim = operations[FeBulk::BulkOps::opJacDet].dim();
        uint result_begin_col = operations[FeBulk::BulkOps::opJacDet].result_col();
        uint input_begin_col = operations[FeBulk::BulkOps::opJac].result_col();
        ArrayDbl &result_vec = op_results(result_begin_col);
        switch (dim) {
            case 1: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 1>> input_mat(op_results.data() + input_begin_col, 3, 1);
                result_vec = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, 1>>(input_mat);
                break;
            }
            case 2: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> input_mat(op_results.data() + input_begin_col, 3, 2);
                result_vec = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, 2>>(input_mat);
                break;
            }
            case 3: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 3>> input_mat(op_results.data() + input_begin_col, 3, 3);
                result_vec = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, 3>>(input_mat);
                break;
            }
        }
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void ptop_jac(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void ptop_jac_det(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
};



namespace FeSide {

/// enum of side operation
enum SideOps
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

} // closing namespace FeSide



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
	static inline void elop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void elop_jac(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        uint dim = operations[FeSide::SideOps::opJac].dim();
        uint result_begin_col = operations[FeSide::SideOps::opJac].result_col();
        uint input_begin_col = operations[FeSide::SideOps::opElCoords].result_col();
        switch (dim) {
            // no evaluation for dim=0, shape of Jacobian (spacedim,0)
            case 1: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 1>> result_mat(op_results.data() + result_begin_col, 3, 1);
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> input_mat(op_results.data() + input_begin_col, 3, 2);
                result_mat = eigen_tools::jacobian<3,1>(input_mat);
                break;
            }
            case 2: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> result_mat(op_results.data() + result_begin_col, 3, 2);
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 3>> input_mat(op_results.data() + input_begin_col, 3, 3);
                result_mat = eigen_tools::jacobian<3,2>(input_mat);
                break;
            }
        }
    }
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results) {
        // result double, input matrix(spacedim, dim)
        uint dim = operations[FeSide::SideOps::opJacDet].dim();
        uint result_begin_col = operations[FeSide::SideOps::opJacDet].result_col();
        uint input_begin_col = operations[FeSide::SideOps::opJac].result_col();
        ArrayDbl &result_vec = op_results(result_begin_col);
        switch (dim) {
            case 0: {
                for (uint i=0;i<300; ++i)
                    result_vec(i) = 1.0;
                break;
            }
            case 1: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 1>> input_mat(op_results.data() + input_begin_col, 3, 1);
                result_vec = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, 1>>(input_mat);
                break;
            }
            case 2: {
                Eigen::Map<Eigen::Matrix<ArrayDbl, 3, 2>> input_mat(op_results.data() + input_begin_col, 3, 2);
                result_vec = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, 2>>(input_mat);
                break;
            }
        }
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void ptop_jac(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
    static inline void ptop_jac_det(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results) {
        // Implement
    }
};



#endif /* PATCH_POINT_VALUES_HH_ */
