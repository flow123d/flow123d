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
#include "fem/element_values.hh"


template<unsigned int spacedim> class PatchFEValues;
template<unsigned int spacedim> class ElOp;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;




/// Type for conciseness
typedef void (*ReinitFunction)(std::vector<ElOp<3>> &, TableDbl &, TableInt &);


namespace FeBulk {
    /// enum of bulk operation
    enum BulkOps
    {
        opElCoords,
        opJac,
        opJacDet,
        opExpdCoords,
        opExpdJac,
        opExpdJacDet,
        opCoords,
        opWeights,
        opJxW
    };
}


namespace FeSide {
    /// enum of side operation
    enum SideOps
    {
        opElCoords,
        opJac,
        opJacDet,
        opExpdCoords,
        opExpdJac,
        opExpdJacDet,
        opCoords,
        opWeights,
        opJxW
    };
}



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
    PatchPointValues(uint dim)
    : dim_(dim), n_columns_(0), elements_map_(300, 0), points_map_(300, 0) {}

    /**
     * Initialize object, set number of columns (quantities) in tables.
     *
     * Number of columns of int_vals_ table is passed by argument, number of columns
     * of other tables is given by n_columns_ value.
     */
    void initialize(Quadrature &quad, uint int_cols) {
        this->reset();

    	point_vals_.resize(n_columns_);
    	int_vals_.resize(int_cols);

    	ref_elm_values_ = std::make_shared<RefElementValues<spacedim> >(quad, dim_);
    	ref_elm_values_->ref_initialize(quad, dim_);
    }

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

    /// Resize data tables
    void resize_tables(uint n_points) {
        eigen_tools::resize_table(point_vals_, n_points);
        eigen_tools::resize_table(int_vals_, n_points);
    }

    /// Register element, add to point_vals_ table
    uint register_element(arma::mat coords, uint element_patch_idx) {
        uint res_column = operations_[FeBulk::BulkOps::opElCoords].result_col();
        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
                point_vals_(res_column)(n_elems_) = coords(i_row, i_col);
                ++res_column;
            }

        elements_map_[element_patch_idx] = n_elems_;
        return n_elems_++;
    }

    /// Register side, add to point_vals_ table
    uint register_side(arma::mat coords) {
        uint res_column = operations_[FeSide::SideOps::opElCoords].result_col();
        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
                point_vals_(res_column)(n_elems_) = coords(i_row, i_col);
                ++res_column;
            }

        return n_elems_++;
    }

    /// Register bulk point, add to int_vals_ table
    uint register_bulk_point(uint elem_table_row, uint value_patch_idx, uint elem_idx) {
        int_vals_(0)(n_points_) = value_patch_idx;
        int_vals_(1)(n_points_) = elem_table_row;
        int_vals_(2)(n_points_) = elem_idx;

        points_map_[value_patch_idx] = n_points_;
        return n_points_++;
    }

    /// Register side point, add to int_vals_ table
    uint register_side_point(uint elem_table_row, uint value_patch_idx, uint elem_idx, uint side_idx) {
        int_vals_(0)(n_points_) = value_patch_idx;
        int_vals_(1)(n_points_) = elem_table_row;
        int_vals_(2)(n_points_) = elem_idx;
        int_vals_(3)(n_points_) = side_idx;

        points_map_[value_patch_idx] = n_points_;
        return n_points_++;
    }

    /// Add accessor to operations_ vector - temporary method
    ElOp<spacedim> &add_accessor(ElOp<spacedim> op_accessor) {
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }

    /// Add accessor to operations_ vector
    ElOp<spacedim> &make_new_op(std::initializer_list<uint> shape, ReinitFunction reinit_func, std::vector<uint> input_ops_vec) {
    	ElOp<spacedim> op_accessor(this->dim_, shape, this->n_columns_, reinit_func, input_ops_vec);
    	this->n_columns_ += op_accessor.n_comp();
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }

    /// Add accessor to operations_ vector
    ElOp<spacedim> &make_expansion(ElOp<spacedim> &el_op, std::initializer_list<uint> shape, ReinitFunction reinit_func) {
        ElOp<spacedim> op_accessor(this->dim_, shape, el_op.result_col(), reinit_func);
        // shape passed from el_op throws:
        // C++ exception with description "std::bad_alloc" thrown in the test body.
        operations_.push_back(op_accessor);
        return operations_[operations_.size()-1];
    }


    /// Reinit data.
    void reinit_patch() {
        std::cout << "Reinit_patch: " << dim_ << std::endl;
        if (n_elems_ == 0) return; // skip if tables are empty
        // precompute data on point_vals_ table
        std::cout << "Reinit_patch A" << std::endl;
        for (uint i=0; i<operations_.size(); ++i)
            operations_[i].reinit_funtionc(operations_, point_vals_, int_vals_);
        std::cout << "Reinit_patch B" << std::endl;

//        // copy data from el_vals_ to point_vals_
//        std::vector<uint> copied; // list of columns that will be copied
//        for (uint i_op=0; i_op<operations_.size(); ++i_op)
//            if (operations_[i_op].copy_vals()) {
//            	uint res_col = operations_[i_op].result_col();
//            	for (uint i_col=res_col; i_col<res_col+operations_[i_op].n_comp(); ++i_col)
//            		copied.push_back(i_col);
//            }
//        for (uint i_pt=0; i_pt<n_points_; ++i_pt) {
//            uint el_table_idx = int_vals_(1)(i_pt);
//            for (uint i_q=0; i_q<copied.size(); ++i_q)
//                point_vals_(copied[i_q])(i_pt) = el_vals_(copied[i_q])(el_table_idx);
//        }
//
//        // precompute data on point_vals_ table
//        for (uint i=0; i<operations_.size(); ++i)
//            operations_[i].reinit_points(operations_, point_vals_);
    }

    /// Temporary development method
    void print(bool points, bool ints) const {
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
     *  1: Row of element/side in point_vals_ table in registration step (before expansion)
     *  2: Element idx in Mesh
     *  3: Index of side in element (column is allocated only for side point table)
     * Number of used rows is given by n_points_.
     */
    TableInt int_vals_;

    /// Vector of all defined operations
    std::vector<ElOp<spacedim>> operations_;

    uint dim_;                        ///< Dimension
    uint n_columns_;                  ///< Number of columns of \p point_vals table
    uint n_points_;                   ///< Number of points in patch
    uint n_elems_;                    ///< Number of elements in patch

    std::vector<uint> elements_map_;  ///< Map of element patch indices to el_vals_ table
    std::vector<uint> points_map_;    ///< Map of point patch indices  to point_vals_ and int_vals_ tables

    /// Auxiliary object for calculation of element-dependent data.
    std::shared_ptr<RefElementValues<spacedim> > ref_elm_values_;

    friend class PatchFEValues<spacedim>;
    friend class ElOp<spacedim>;
};


/**
 * @brief Class represents FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
    /// Constructor
    ElOp(uint dim, std::initializer_list<uint> shape, uint result_col, ReinitFunction reinit_fn = nullptr,
            std::vector<uint> input_ops = {})
    : dim_(dim), shape_(shape), result_col_(result_col), input_ops_(input_ops), reinit_func(reinit_fn)
    {}

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

    /// Getter of input_ops_
    inline const std::vector<uint> &input_ops() const {
        return input_ops_;
    }

    /// Call reinit function on element table if function is defined
    inline void reinit_funtionc(std::vector<ElOp<spacedim>> &operations, TableDbl &data_table, TableInt &int_table) {
    	if (reinit_func != nullptr) reinit_func(operations, data_table, int_table);
    }

//    inline Scalar scalar_val(uint point_idx) const {
//        return point_vals_->point_vals_[result_col_][point_idx];
//    }
//
//    inline Vector vector_val(uint point_idx) const {
//        Vector val;
//        for (uint i=0; i<3; ++i)
//            val(i) = point_vals_->point_vals_[result_col_+i][point_idx];
//        return val;
//    }
//
//    inline Tensor tensor_val(uint point_idx) const {
//        Tensor val;
//        for (uint i=0; i<3; ++i)
//            for (uint j=0; j<3; ++j)
//                val(i,j) = point_vals_->point_vals_[result_col_+3*i+j][point_idx];
//        return val;
//    }


protected:
    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    uint result_col_;                         ///< First column to scalar, vector or matrix result
    std::vector<uint> input_ops_;             ///< Indices of operations in PatchPointValues::operations_ vector on which ElOp is depended

    ReinitFunction reinit_func;               ///< Pointer to patch reinit function of element data table specialized by operation
};


/// Defines common functionality of reinit operations.
struct common_reinit {
	// empty base operation
	static inline void op_base(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // empty
    }

	// expansion operation
    static inline void expand_data(ElOp<3> &op, TableDbl &op_results, TableInt &el_table) {
        uint row_begin = op.result_col();
        uint row_end = row_begin + op.n_comp();
        uint size = op_results(row_begin).rows();
        for (int i_pt=size-1; i_pt>=0; --i_pt) {
            uint el_table_idx = el_table(1)(i_pt);
            for (uint i_q=row_begin; i_q<row_end; ++i_q)
                op_results(i_q)(i_pt) = op_results(i_q)(el_table_idx);
        }
    }
};

/// Defines reinit operations on bulk points.
struct bulk_reinit {
	// element operations
    static inline void elop_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeBulk::BulkOps::opJac];
        uint result_begin_col = op.result_col();
        uint input_begin_col = operations[op.input_ops()[0]].result_col();
        switch (op.dim()) {
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
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeBulk::BulkOps::opJacDet];
        uint result_begin_col = op.result_col();
        uint input_begin_col = operations[op.input_ops()[0]].result_col();
        ArrayDbl &result_vec = op_results(result_begin_col);
        switch (op.dim()) {
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

    // expansion operations
    static inline void expd_coords(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpdCoords];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpdJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpdJacDet];
        common_reinit::expand_data(op, op_results, el_table);
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_weights(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_JxW(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
};



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
    static inline void elop_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeSide::SideOps::opJac];
        uint result_begin_col = op.result_col();
        uint input_begin_col = operations[op.input_ops()[0]].result_col();
        switch (op.dim()) {
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
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeSide::SideOps::opJacDet];
        uint result_begin_col = op.result_col();
        uint input_begin_col = operations[op.input_ops()[0]].result_col();
        ArrayDbl &result_vec = op_results(result_begin_col);
        switch (op.dim()) {
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

    // expansion operations
    static inline void expd_coords(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpdCoords];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpdJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpdJacDet];
        common_reinit::expand_data(op, op_results, el_table);
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_weights(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_JxW(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
};


namespace FeBulk {

    /// Bulk data specialization, order of item in operations_ vector corresponds to the BulkOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        /// Constructor fills operations_ vector.
        PatchPointValues(uint dim)
        : ::PatchPointValues<spacedim>(dim) {
            // First step: adds element values operations
            auto &el_coords = this->make_new_op(
                    {spacedim, this->dim_+1}, &common_reinit::op_base, {} );

            auto &el_jac = this->make_new_op(
                    {spacedim, this->dim_}, &bulk_reinit::elop_jac, {BulkOps::opElCoords} );

            auto &el_jac_det = this->make_new_op(
                    {1}, &bulk_reinit::elop_jac_det, {BulkOps::opJac} );

            // Second step: adds expand operations (element values to point values)
            this->make_expansion(
                    el_coords, {spacedim, this->dim_+1}, &bulk_reinit::expd_coords );

            this->make_expansion(
                    el_jac, {spacedim, this->dim_}, &bulk_reinit::expd_jac );

            this->make_expansion(
                    el_jac_det, {1}, &bulk_reinit::expd_jac_det );

            // Third step: adds point values operations
            ElOp<spacedim> &coords = this->add_accessor(
            		ElOp<spacedim>(this->dim_, {spacedim}, this->n_columns_,
            				&bulk_reinit::ptop_coords) );
            this->n_columns_ += coords.n_comp();

            ElOp<spacedim> &weights = this->add_accessor(
            		ElOp<spacedim>(this->dim_, {1}, this->n_columns_,
            				&bulk_reinit::ptop_weights) );
            this->n_columns_ += weights.n_comp();

            ElOp<spacedim> &JxW = this->add_accessor(
            		ElOp<spacedim>(this->dim_, {1}, this->n_columns_,
            				&bulk_reinit::ptop_JxW, {BulkOps::opWeights}) );
            this->n_columns_ += JxW.n_comp();
            std::cout << "Constructor B" << std::endl;
        }
    };

} // closing namespace FeBulk



namespace FeSide {

/// Bulk Side specialization, order of item in operations_ vector corresponds to the SideOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        /// Constructor fills operations_ vector.
        PatchPointValues(uint dim)
        : ::PatchPointValues<spacedim>(dim) {
            // First step: adds element values operations
            auto &el_coords = this->make_new_op( {spacedim, this->dim_+1}, &common_reinit::op_base, {} );

            auto &el_jac = this->make_new_op( {spacedim, this->dim_}, &side_reinit::elop_jac, {SideOps::opElCoords} );

            auto &el_jac_det = this->make_new_op( {1}, &side_reinit::elop_jac_det, {SideOps::opJac} );

            // Second step: adds expand operations (element values to point values)
            this->make_expansion(
                    el_coords, {spacedim, this->dim_+1}, &side_reinit::expd_coords );

            this->make_expansion(
                    el_jac, {spacedim, this->dim_}, &side_reinit::expd_jac );

            this->make_expansion(
                    el_jac_det, {1}, &side_reinit::expd_jac_det );

            // Third step: adds point values operations
            ElOp<spacedim> &coords = this->add_accessor( ElOp<spacedim>(this->dim_, {spacedim}, this->n_columns_, &side_reinit::ptop_coords) );
            this->n_columns_ += coords.n_comp();

            ElOp<spacedim> &weights = this->add_accessor( ElOp<spacedim>(this->dim_, {1}, this->n_columns_, &side_reinit::ptop_weights) );
            this->n_columns_ += weights.n_comp();

            ElOp<spacedim> &JxW = this->add_accessor( ElOp<spacedim>(this->dim_, {1}, this->n_columns_, &side_reinit::ptop_JxW, {SideOps::opWeights}) );
            this->n_columns_ += JxW.n_comp();
        }
    };

} // closing namespace FeSide



#endif /* PATCH_POINT_VALUES_HH_ */
