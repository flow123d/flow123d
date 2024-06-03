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
#include "quadrature/quadrature_lib.hh"


template<unsigned int spacedim> class PatchFEValues;
template<unsigned int spacedim> class ElOp;
template<unsigned int dim> class BulkValues;
template<unsigned int dim> class SideValues;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;




/// Type for conciseness
using ReinitFunction = std::function<void(std::vector<ElOp<3>> &, TableDbl &, TableInt &)>;


namespace FeBulk {
    /// enum of bulk operation
    enum BulkOps
    {
        /// operations evaluated on elements
        opElCoords,         ///< coordinations of all nodes of element
        opJac,              ///< Jacobian of element
        opInvJac,           ///< inverse Jacobian
        opJacDet,           ///< determinant of Jacobian
		/// operation executed expansion to quadrature point (value of element > values on quadrature points)
        opExpansionCoords,  ///< expands coordinates
        opExpansionJac,     ///< expands Jacobian
        opExpansionInvJac,  ///< expands inverse Jacobian
        opExpansionJacDet,  ///< expands Jacobian determinant
        /// operations evaluated on quadrature points
        opCoords,           ///< coordinations of quadrature point
        opWeights,          ///< weight of quadrature point
        opJxW               ///< JxW value of quadrature point
    };
}


namespace FeSide {
    /// enum of side operation
    enum SideOps
    {
        /// operations evaluated on elements
        opElCoords,             ///< coordinations of all nodes of element
        opElJac,                ///< Jacobian of element
		opElInvJac,             ///< inverse Jacobian of element
        /// operations evaluated on sides
        opSideCoords,           ///< coordinations of all nodes of side
        opSideJac,              ///< Jacobian of element
        opSideJacDet,           ///< determinant of Jacobian of side
		/// operation executed expansion to quadrature point (value of element / side > values on quadrature points)
        opExpansionElCoords,    ///< expands coordinates on element
        opExpansionElJac,       ///< expands Jacobian on element
        opExpansionElInvJac,    ///< expands inverse Jacobian on element
        opExpansionSideCoords,  ///< expands coordinates on side
        opExpansionSideJac,     ///< expands Jacobian on side
        opExpansionSideJacDet,  ///< expands Jacobian determinant on side
        /// operations evaluated on quadrature points
        opCoords,               ///< coordinations of quadrature point
        opWeights,              ///< weight of quadrature point
        opJxW,                  ///< JxW value of quadrature point
		opNormalVec             ///< normal vector of quadrature point
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
    : dim_(dim), n_rows_(0), elements_map_(300, 0), points_map_(300, 0) {}

    /**
     * Initialize object, set number of columns (quantities) in tables.
     *
     * Number of columns of int_vals_ table is passed by argument, number of columns
     * of other tables is given by n_rows_ value.
     */
    void initialize(uint int_cols) {
        this->reset();

    	point_vals_.resize(n_rows_);
    	int_vals_.resize(int_cols);
    }

    /// Reset number of rows (points and elements)
    inline void reset() {
        n_points_ = 0;
        n_elems_ = 0;
    }

    /// Getter of dim_
    inline uint dim() const {
        return dim_;
    }

    /// Getter of n_rows_
    inline uint n_rows() const {
        return n_rows_;
    }

    /// Getter of n_elems_
    inline uint n_elems() const {
        return n_elems_;
    }

    /// Getter of n_points_
    inline uint n_points() const {
        return n_points_;
    }

    /// Return quadrature
    Quadrature *get_quadrature() const {
        return quad_;
    }

    /// Resize data tables
    void resize_tables(uint n_points) {
        eigen_tools::resize_table(point_vals_, n_points);
        eigen_tools::resize_table(int_vals_, n_points);
    }

    /// Register element, add to point_vals_ table
    uint register_element(arma::mat coords, uint element_patch_idx) {
        uint res_column = operations_[FeBulk::BulkOps::opElCoords].result_row();
        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
                point_vals_(res_column)(n_elems_) = coords(i_row, i_col);
                ++res_column;
            }

        elements_map_[element_patch_idx] = n_elems_;
        return n_elems_++;
    }

    /// Register side, add to point_vals_ table
    uint register_side(arma::mat elm_coords, arma::mat side_coords) {
        uint res_column = operations_[FeSide::SideOps::opElCoords].result_row();
        for (uint i_col=0; i_col<elm_coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<elm_coords.n_rows; ++i_row) {
                point_vals_(res_column)(n_elems_) = elm_coords(i_row, i_col);
                ++res_column;
            }

        res_column = operations_[FeSide::SideOps::opSideCoords].result_row();
        for (uint i_col=0; i_col<side_coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<side_coords.n_rows; ++i_row) {
                point_vals_(res_column)(n_elems_) = side_coords(i_row, i_col);
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

    /// Add accessor to operations_ vector
    ElOp<spacedim> &make_new_op(std::initializer_list<uint> shape, ReinitFunction reinit_f, std::vector<uint> input_ops_vec) {
    	ElOp<spacedim> op_accessor(this->dim_, shape, this->n_rows_, reinit_f, input_ops_vec);
    	this->n_rows_ += op_accessor.n_comp();
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }

    /// Add accessor to operations_ vector
    ElOp<spacedim> &make_expansion(ElOp<spacedim> &el_op, std::initializer_list<uint> shape, ReinitFunction reinit_f) {
        ElOp<spacedim> op_accessor(this->dim_, shape, el_op.result_row(), reinit_f);
        // shape passed from el_op throws:
        // C++ exception with description "std::bad_alloc" thrown in the test body.
        operations_.push_back(op_accessor);
        return operations_[operations_.size()-1];
    }

    /// Add accessor to operations_ vector
    ElOp<spacedim> &make_fe_op(std::initializer_list<uint> shape, ReinitFunction reinit_f, std::vector<uint> input_ops_vec, uint n_dofs) {
    	ElOp<spacedim> op_accessor(this->dim_, shape, this->n_rows_, reinit_f, input_ops_vec, n_dofs);
    	this->n_rows_ += op_accessor.n_comp() * n_dofs;
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }


    /// Reinit data.
    void reinit_patch() {
        if (n_elems_ == 0) return; // skip if tables are empty
        // Reinit patch data of all operation
        for (uint i=0; i<operations_.size(); ++i)
            operations_[i].reinit_function(operations_, point_vals_, int_vals_);
    }

    inline Scalar scalar_val(uint result_row, uint point_idx) const {
        return point_vals_(result_row)(points_map_[point_idx]);
    }

    inline Vector vector_val(uint result_row, uint point_idx) const {
        Vector val;
        for (uint i=0; i<3; ++i)
            val(i) = point_vals_(result_row+i)(points_map_[point_idx]);
        return val;
    }

    inline Tensor tensor_val(uint result_row, uint point_idx) const {
        Tensor val;
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = point_vals_(result_row+3*i+j)(points_map_[point_idx]);
        return val;
    }

    /// Temporary development method
    void print_data_tables(ostream& stream, bool points, bool ints) const {
        if (points) {
            stream << "Point vals: " << point_vals_.rows() << " - " << point_vals_.cols() << std::endl;
	        for (uint i_row=0; i_row<n_points_; ++i_row) {
                for (uint i_col=0; i_col<n_rows_; ++i_col)
                	stream << point_vals_(i_col)(i_row) << " ";
                stream << std::endl;
            }
            stream << std::endl;
        }
        if (ints) {
            stream << "Int vals: " << int_vals_.rows() << " - " << int_vals_.cols() << std::endl;
	        for (uint i_row=0; i_row<n_points_; ++i_row) {
                for (uint i_col=0; i_col<3; ++i_col)
                	stream << int_vals_(i_col)(i_row) << " ";
                stream << std::endl;
            }
            stream << std::endl;
        }
    }

    /// Temporary development method
    void print_operations(ostream& stream, uint bulk_side) const {
        std::vector< std::vector<std::string> > op_names =
        {
            { "el_coords", "jacobian", "inv_jac", "jac_det", "exp_coords", "exp_jacobian", "exp_in_jac", "exp_jac_det", "pt_coords", "weights",
              "JxW", "", "", "", "", "" },
            { "el_coords", "el_jac", "el_inv_jac", "side_coords", "side_jac", "side_jac_det", "exp_el_coords", "exp_el_jac", "exp_el_inv_jac",
              "exp_side_coords", "exp_side_jac", "exp_side_jac_det", "pt_coords", "weights", "JxW", "normal_vec", "", "", "", "", "" }
        };
        stream << std::setfill(' ') << " Operation" << setw(12) << "" << "Shape" << setw(2) << ""
                << "Result row" << setw(2) << "" << "n DOFs" << setw(2) << "" << "Input operations" << endl;
        for (uint i=0; i<operations_.size(); ++i) {
            stream << " " << std::left << setw(20) << op_names[bulk_side][i] << "" << " " << setw(6) << operations_[i].format_shape() << "" << " "
                << setw(11) << operations_[i].result_row() << "" << " " << setw(7) << operations_[i].n_dofs() << "" << " ";
            auto &input_ops = operations_[i].input_ops();
            for (auto i_o : input_ops) stream << op_names[bulk_side][i_o] << " ";
            stream << std::endl;
        }
    }

protected:
    /**
     * Store data of bulk or side quadrature points of one dimension
     *
     * Number of columns is given by n_rows_, number of used rows by n_points_.
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
    uint n_rows_;                     ///< Number of columns of \p point_vals table
    uint n_points_;                   ///< Number of points in patch
    uint n_elems_;                    ///< Number of elements in patch
    Quadrature *quad_;                ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;  ///< Map of element patch indices to el_vals_ table
    std::vector<uint> points_map_;    ///< Map of point patch indices  to point_vals_ and int_vals_ tables

    friend class PatchFEValues<spacedim>;
    friend class ElOp<spacedim>;
    template<unsigned int dim>
    friend class BulkValues;
    template<unsigned int dim>
    friend class SideValues;
    template<unsigned int dim>
    friend class JoinValues;
};

/**
 * @brief Class represents FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
    /// Constructor
    ElOp(uint dim, std::initializer_list<uint> shape, uint result_row, ReinitFunction reinit_f, std::vector<uint> input_ops = {}, uint n_dofs = 1)
    : dim_(dim), shape_(shape), result_row_(result_row), input_ops_(input_ops), n_dofs_(n_dofs), reinit_func(reinit_f)
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

    /// Getter of result_row_
    inline uint result_row() const {
        return result_row_;
    }

    /// Getter of n_dofs_
    inline uint n_dofs() const {
        return n_dofs_;
    }

    /// Getter of input_ops_
    inline const std::vector<uint> &input_ops() const {
        return input_ops_;
    }

    /// Getter of shape_
    inline const std::vector<uint> &shape() const {
        return shape_;
    }

    /// Format shape to string
    inline std::string format_shape() const {
        stringstream ss;
        ss << shape_[0];
        if (shape_.size() > 1) ss << "x" << shape_[1];
        return ss.str();
    }

    /// Call reinit function on element table if function is defined
    inline void reinit_function(std::vector<ElOp<spacedim>> &operations, TableDbl &data_table, TableInt &int_table) {
        reinit_func(operations, data_table, int_table);
    }

    /// Return map referenced Eigen::Matrix of given dimension
    template<unsigned int dim1, unsigned int dim2>
    Eigen::Map<Eigen::Matrix<ArrayDbl, dim1, dim2>> value(TableDbl &op_results, uint i_dof = 0) const {
        return Eigen::Map<Eigen::Matrix<ArrayDbl, dim1, dim2>>(op_results.data() + result_row_ + i_dof * n_comp(), dim1, dim2);
    }


protected:
    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    uint result_row_;                         ///< First row to scalar, vector or matrix result
    std::vector<uint> input_ops_;             ///< Indices of operations in PatchPointValues::operations_ vector on which ElOp is depended
    uint n_dofs_;                             ///< Number of DOFs of FE operations (or 1 in case of element operations)

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
        uint row_begin = op.result_row();
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
    template<unsigned int dim>
    static inline void elop_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeBulk::BulkOps::opJac];
        auto jac_value = op.value<3, dim>(op_results);
        auto coords_value = operations[ op.input_ops()[0] ].value<3, dim+1>(op_results);
        jac_value = eigen_tools::jacobian<3,dim>(coords_value);
    }
    template<unsigned int dim>
    static inline void elop_inv_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeBulk::BulkOps::opInvJac];
        auto inv_jac_value = op.value<dim, 3>(op_results);
        auto jac_value = operations[ op.input_ops()[0] ].value<3, dim>(op_results);
        inv_jac_value = eigen_tools::inverse<3, dim>(jac_value);
    }
    template<unsigned int dim>
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeBulk::BulkOps::opJacDet];
        auto &jac_det_value = op_results(op.result_row());
        auto jac_value = operations[ op.input_ops()[0] ].value<3, dim>(op_results);
        jac_det_value = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, dim>>(jac_value).array().abs();
    }

    // expansion operations
    static inline void expd_coords(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpansionCoords];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpansionJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_inv_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpansionInvJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opExpansionJacDet];
        common_reinit::expand_data(op, op_results, el_table);
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_weights(std::vector<ElOp<3>> &operations, TableDbl &op_results, std::vector<double> point_weights) {
        auto &op = operations[FeBulk::BulkOps::opWeights];
        ArrayDbl &result_row = op_results( op.result_row() );
        auto n_points = point_weights.size();
        for (uint i=0; i<result_row.rows(); ++i)
            result_row(i) = point_weights[i%n_points];
    }
    static inline void ptop_JxW(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        auto &op = operations[FeBulk::BulkOps::opJxW];
        ArrayDbl &weights_row = op_results( operations[op.input_ops()[0]].result_row() );
        ArrayDbl &jac_det_row = op_results( operations[op.input_ops()[1]].result_row() );
        ArrayDbl &result_row = op_results( op.result_row() );
        result_row = jac_det_row * weights_row;
    }
    static inline void ptop_scalar_shape(std::vector<ElOp<3>> &operations, TableDbl &op_results,
            std::vector< std::vector<double> > shape_values, uint scalar_shape_op_idx) {
        auto &op = operations[scalar_shape_op_idx];
        uint n_points = shape_values.size();

        for (uint i_row=0; i_row<shape_values[0].size(); ++i_row) {
            ArrayDbl &result_row = op_results( op.result_row()+i_row );
            for (uint i_pt=0; i_pt<result_row.rows(); ++i_pt)
                result_row(i_pt) = shape_values[i_pt % n_points][i_row];
        }
    }
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(std::vector<ElOp<3>> &operations, TableDbl &op_results,
            std::vector< std::vector<arma::mat> > ref_shape_grads, uint scalar_shape_grads_op_idx) {
        auto &op = operations[scalar_shape_grads_op_idx];
        uint n_points = ref_shape_grads.size();
        uint n_dofs = ref_shape_grads[0].size();

        Eigen::Vector<ArrayDbl, Eigen::Dynamic> ref_shape_grads_expd;
        ref_shape_grads_expd.resize(dim*n_dofs);
        for (uint i=0; i<ref_shape_grads_expd.rows(); ++i)
        	ref_shape_grads_expd(i).resize(op_results(0).rows());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_c=0; i_c<dim; ++i_c) {
                ArrayDbl &shape_grad_row = ref_shape_grads_expd(i_dof*dim+i_c);
                for (uint i_pt=0; i_pt<shape_grad_row.rows(); ++i_pt)
                    shape_grad_row(i_pt) = ref_shape_grads[i_pt % n_points][i_dof][i_c];
            }

        auto inv_jac_value = operations[ op.input_ops()[0] ].value<dim, 3>(op_results);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            auto shape_grad_value = op.value<3, 1>(op_results, i_dof);
            Eigen::Map<Eigen::Matrix<ArrayDbl, dim, 1>> ref_shape_grads_dof_value(ref_shape_grads_expd.data() + dim*i_dof, dim, 1);
            shape_grad_value = inv_jac_value.transpose() * ref_shape_grads_dof_value;
        }
    }
};



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_el_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeSide::SideOps::opElJac];
        auto jac_value = op.value<3, dim>(op_results);
        auto coords_value = operations[ op.input_ops()[0] ].value<3, dim+1>(op_results);
        jac_value = eigen_tools::jacobian<3,dim>(coords_value);
    }
    template<unsigned int dim>
    static inline void elop_el_inv_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opElInvJac];
        auto inv_jac_value = op.value<dim, 3>(op_results);
        auto jac_value = operations[ op.input_ops()[0] ].value<3, dim>(op_results);
        inv_jac_value = eigen_tools::inverse<3,dim>(jac_value);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeSide::SideOps::opSideJac];
        auto jac_value = op.value<3, dim-1>(op_results);
        auto coords_value = operations[ op.input_ops()[0] ].value<3, dim>(op_results);
        jac_value = eigen_tools::jacobian<3, dim-1>(coords_value);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeSide::SideOps::opSideJacDet];
        ArrayDbl &det_value = op_results( op.result_row() );
        auto jac_value = operations[ op.input_ops()[0] ].value<3, dim-1>(op_results);
        det_value = eigen_tools::determinant<Eigen::Matrix<ArrayDbl, 3, dim-1>>(jac_value).array().abs();
    }

    // expansion operations
    static inline void expd_el_coords(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionElCoords];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_el_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionElJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_el_inv_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionElInvJac];
        common_reinit::expand_data(op, op_results, el_table);
    }

    static inline void expd_sd_coords(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionSideCoords];
        common_reinit::expand_data(op, op_results, el_table);
    }
    template<unsigned int dim>
    static inline void expd_sd_jac(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionSideJac];
        common_reinit::expand_data(op, op_results, el_table);
    }
    static inline void expd_sd_jac_det(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opExpansionSideJacDet];
        common_reinit::expand_data(op, op_results, el_table);
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        // Implement
    }
    static inline void ptop_weights(std::vector<ElOp<3>> &operations, TableDbl &op_results, std::vector<double> point_weights) {
        auto &op = operations[FeSide::SideOps::opWeights];
        ArrayDbl &result_row = op_results( op.result_row() );
        auto n_points = point_weights.size();
        for (uint i=0; i<result_row.rows(); ++i)
            result_row(i) = point_weights[i%n_points];
    }
    static inline void ptop_JxW(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opJxW];
        ArrayDbl &weights_row = op_results( operations[op.input_ops()[0]].result_row() );
        ArrayDbl &jac_det_row = op_results( operations[op.input_ops()[1]].result_row() );
        ArrayDbl &result_row = op_results( op.result_row() );
        result_row = jac_det_row * weights_row;
    }
    template<unsigned int dim>
    static inline void ptop_normal_vec(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table) {
        auto &op = operations[FeSide::SideOps::opNormalVec];
        auto normal_value = op.value<3, 1>(op_results);
        auto inv_jac_mat_value = operations[ op.input_ops()[0] ].value<dim, 3>(op_results);
        normal_value = inv_jac_mat_value.transpose() * RefElement<dim>::normal_vector_array( el_table(3) );

        ArrayDbl norm_vec;
        norm_vec.resize(normal_value(0).rows());
        Eigen::VectorXd A(3);

        for (uint i=0; i<normal_value(0).rows(); ++i) {
            A(0) = normal_value(0)(i);
            A(1) = normal_value(1)(i);
            A(2) = normal_value(2)(i);
            norm_vec(i) = A.norm();
        }
        for (uint i=0; i<3; ++i) {
        	normal_value(i) /= norm_vec;
        }
    }
    static inline void ptop_scalar_shape(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table,
    		std::vector< std::vector< std::vector<double> > > shape_values, uint scalar_shape_op_idx) {
        auto &op = operations[scalar_shape_op_idx];
        uint n_points = shape_values[0].size();

        for (uint i_row=0; i_row<shape_values[0][0].size(); ++i_row) {
            ArrayDbl &result_row = op_results( op.result_row()+i_row );
            for (uint i_pt=0; i_pt<result_row.rows(); ++i_pt)
                result_row(i_pt) = shape_values[el_table(3)(i_pt)][i_pt % n_points][i_row];
        }
    }
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(std::vector<ElOp<3>> &operations, TableDbl &op_results, TableInt &el_table,
            std::vector< std::vector< std::vector<arma::mat> > > ref_shape_grads, uint scalar_shape_grads_op_idx) {
        auto &op = operations[scalar_shape_grads_op_idx];
        uint n_points = ref_shape_grads[0].size();
        uint n_dofs = ref_shape_grads[0][0].size();

        Eigen::Vector<ArrayDbl, Eigen::Dynamic> ref_shape_grads_expd;
        ref_shape_grads_expd.resize(dim*n_dofs);
        for (uint i=0; i<ref_shape_grads_expd.rows(); ++i)
        	ref_shape_grads_expd(i).resize(op_results(0).rows());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_c=0; i_c<dim; ++i_c) {
                ArrayDbl &shape_grad_row = ref_shape_grads_expd(i_dof*dim+i_c);
                for (uint i_pt=0; i_pt<shape_grad_row.rows(); ++i_pt)
                    shape_grad_row(i_pt) = ref_shape_grads[el_table(3)(i_pt)][i_pt % n_points][i_dof][i_c];
            }

        auto inv_jac_value = operations[ op.input_ops()[0] ].value<dim, 3>(op_results);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            auto shape_grad_value = op.value<3, 1>(op_results, i_dof);
            Eigen::Map<Eigen::Matrix<ArrayDbl, dim, 1>> ref_shape_grads_dof_value(ref_shape_grads_expd.data() + dim*i_dof, dim, 1);
            shape_grad_value = inv_jac_value.transpose() * ref_shape_grads_dof_value;
        }
    }
};


// template specialization
template<>
inline void side_reinit::elop_sd_jac<1>(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
}

template<>
inline void side_reinit::elop_sd_jac_det<1>(std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
    auto &op = operations[FeSide::SideOps::opSideJacDet];
    ArrayDbl &result_vec = op_results( op.result_row() );
    for (uint i=0;i<result_vec.size(); ++i) {
        result_vec(i) = 1.0;
    }
}

template<>
inline void side_reinit::expd_sd_jac<1>(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
}



namespace FeBulk {

    /// Bulk data specialization, order of item in operations_ vector corresponds to the BulkOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        /// Constructor
        PatchPointValues(uint dim, uint quad_order)
        : ::PatchPointValues<spacedim>(dim) {
            this->quad_ = new QGauss(dim, 2*quad_order);
            switch (dim) {
            case 1:
            	init<1>();
            	break;
            case 2:
            	init<2>();
            	break;
            case 3:
            	init<3>();
            	break;
            }
        }

    private:
        /// Initialize operations vector
        template<unsigned int dim>
        void init() {
            // First step: adds element values operations
            auto &el_coords = this->make_new_op( {spacedim, this->dim_+1}, &common_reinit::op_base, {} );

            auto &el_jac = this->make_new_op( {spacedim, this->dim_}, &bulk_reinit::elop_jac<dim>, {BulkOps::opElCoords} );

            auto &el_inv_jac = this->make_new_op( {this->dim_, spacedim}, &bulk_reinit::elop_inv_jac<dim>, {BulkOps::opJac} );

            auto &el_jac_det = this->make_new_op( {1}, &bulk_reinit::elop_jac_det<dim>, {BulkOps::opJac} );

            // Second step: adds expand operations (element values to point values)
            this->make_expansion( el_coords, {spacedim, this->dim_+1}, &bulk_reinit::expd_coords );

            this->make_expansion( el_jac, {spacedim, this->dim_}, &bulk_reinit::expd_jac );

            this->make_expansion( el_inv_jac, {this->dim_, spacedim}, &bulk_reinit::expd_inv_jac );

            this->make_expansion( el_jac_det, {1}, &bulk_reinit::expd_jac_det );

            // Third step: adds point values operations
            /*auto &pt_coords =*/ this->make_new_op( {spacedim}, &bulk_reinit::ptop_coords, {} );

            // use lambda reinit function
            std::vector<double> point_weights = this->quad_->get_weights();
            auto lambda_weights = [point_weights](std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
                    bulk_reinit::ptop_weights(operations, op_results, point_weights);
                };
            /*auto &weights =*/ this->make_new_op( {1}, lambda_weights, {} );

            /*auto &JxW =*/ this->make_new_op( {1}, &bulk_reinit::ptop_JxW, {BulkOps::opWeights, BulkOps::opJacDet} );
        }
    };

} // closing namespace FeBulk



namespace FeSide {

/// Bulk Side specialization, order of item in operations_ vector corresponds to the SideOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        /// Constructor
        PatchPointValues(uint dim, uint quad_order)
        : ::PatchPointValues<spacedim>(dim) {
            this->quad_ = new QGauss(dim-1, 2*quad_order);
            switch (dim) {
            case 1:
            	init<1>();
            	break;
            case 2:
            	init<2>();
            	break;
            case 3:
            	init<3>();
            	break;
            }
        }

    private:
        /// Initialize operations vector
        template<unsigned int dim>
        void init() {
			// First step: adds element values operations
            auto &el_coords = this->make_new_op( {spacedim, this->dim_+1}, &common_reinit::op_base, {} );

            auto &el_jac = this->make_new_op( {spacedim, this->dim_}, &side_reinit::elop_el_jac<dim>, {SideOps::opElCoords} );

            auto &el_inv_jac = this->make_new_op( {this->dim_, spacedim}, &side_reinit::elop_el_inv_jac<dim>, {SideOps::opElJac} );

            auto &sd_coords = this->make_new_op( {spacedim, this->dim_}, &common_reinit::op_base, {} );

            auto &sd_jac = this->make_new_op( {spacedim, this->dim_-1}, &side_reinit::elop_sd_jac<dim>, {SideOps::opSideCoords} );

            auto &sd_jac_det = this->make_new_op( {1}, &side_reinit::elop_sd_jac_det<dim>, {SideOps::opSideJac} );

            // Second step: adds expand operations (element values to point values)
            this->make_expansion( el_coords, {spacedim, this->dim_+1}, &side_reinit::expd_el_coords );

            this->make_expansion( el_jac, {spacedim, this->dim_}, &side_reinit::expd_el_jac );

            this->make_expansion( el_inv_jac, {this->dim_, spacedim}, &side_reinit::expd_el_inv_jac );

            this->make_expansion( sd_coords, {spacedim, this->dim_}, &side_reinit::expd_sd_coords );

            this->make_expansion( sd_jac, {spacedim, this->dim_-1}, &side_reinit::expd_sd_jac<dim> );

            this->make_expansion( sd_jac_det, {1}, &side_reinit::expd_sd_jac_det );

            // Third step: adds point values operations
            /*auto &coords =*/ this->make_new_op( {spacedim}, &side_reinit::ptop_coords, {} );

            // use lambda reinit function
            std::vector<double> point_weights = this->quad_->get_weights();
            auto lambda_weights = [point_weights](std::vector<ElOp<3>> &operations, TableDbl &op_results, FMT_UNUSED TableInt &el_table) {
                    side_reinit::ptop_weights(operations, op_results, point_weights);
                };
            /*auto &weights =*/ this->make_new_op( {1}, lambda_weights, {} );

            /*auto &JxW =*/ this->make_new_op( {1}, &side_reinit::ptop_JxW, {SideOps::opWeights, SideOps::opSideJacDet} );

            /*auto &normal_vec =*/ this->make_new_op( {spacedim}, &side_reinit::ptop_normal_vec<dim>, {SideOps::opElInvJac} );
        }
    };

} // closing namespace FeSide



#endif /* PATCH_POINT_VALUES_HH_ */
