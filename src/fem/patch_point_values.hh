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
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"


template<unsigned int spacedim> class PatchFEValues;
template<unsigned int spacedim> class ElOp;
template <class ValueType> class ElQ;
template<unsigned int dim> class BulkValues;
template<unsigned int dim> class SideValues;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;




/// Type for conciseness
using ReinitFunction = std::function<void(std::vector<ElOp<3>> &, IntTableArena &)>;


namespace FeBulk {
    /**
     * Enumeration of element bulk operations
     *
     * Operations are stored in fix order. Order in enum is equal to order
     * in PatchPointVale::operations_ vector. FE operations are added dynamically
     * by request of user.
     */
    enum BulkOps
    {
        /// operations evaluated on elements
        opElCoords,         ///< coordinations of all nodes of element
        opJac,              ///< Jacobian of element
        opInvJac,           ///< inverse Jacobian
        opJacDet,           ///< determinant of Jacobian
        /// operations evaluated on quadrature points
        opCoords,           ///< coordinations of quadrature point
        opWeights,          ///< weight of quadrature point
        opJxW               ///< JxW value of quadrature point
    };
}


namespace FeSide {
    /**
     * Enumeration of element side operations
     *
     * Operations are stored in fix order. Order in enum is equal to order
     * in PatchPointVale::operations_ vector. FE operations are added dynamically
     * by request of user.
     */
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
        /// operations evaluated on quadrature points
        opCoords,               ///< coordinations of quadrature point
        opWeights,              ///< weight of quadrature point
        opJxW,                  ///< JxW value of quadrature point
		opNormalVec             ///< normal vector of quadrature point
    };
}


/**
 * Enumeration of FE operations.
 *
 * PatchFEValues::fe_ops_indices_ vector holds indices of FE operations
 * in operations_ vector or value of invalid operation if operations is
 * not defined.
 */
enum FEOps
{
    opScalarShape,         ///< Scalar shape operation
    opVectorShape,         ///< Vector shape operation
    opGradScalarShape,     ///< Scalar shape gradient
    opGradVectorShape,     ///< Vector shape gradient
    opVectorSymGrad,       ///< Vector symmetric gradient
    opVectorDivergence,    ///< Vector divergence
    OpNItems               ///< Holds number of valid FE operations and value of invalid FE operation
};


/// Distinguishes operations by type and size of output rows
enum OpSizeType
{
	elemOp,      ///< operation is evaluated on elements or sides
	pointOp,     ///< operation is evaluated on quadrature points
	fixedSizeOp  ///< operation has fixed size and it is filled during initialization
};



/**
 * v Class for storing FE data of quadrature points on one patch.
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
     * @param dim Set dimension
     */
    PatchPointValues(uint dim, AssemblyArena &asm_arena)
    : fe_ops_indices_(FEOps::OpNItems, FEOps::OpNItems), dim_(dim), elements_map_(300, 0),
      points_map_(300, 0), asm_arena_(asm_arena) {
        reset();
    }

    /**
     * Initialize object, set number of columns (quantities) in tables.
     */
    void initialize() {
        this->reset();
        int_table_.resize(int_sizes_.size());
    }

    inline void init_finalize(PatchArena *patch_arena) {
        patch_arena_ = patch_arena;
    }

    /// Reset number of columns (points and elements)
    inline void reset() {
        n_points_ = 0;
        n_elems_ = 0;
        i_elem_ = 0;
    }

    /// Getter for dim_
    inline uint dim() const {
        return dim_;
    }

    /// Getter for n_elems_
    inline uint n_elems() const {
        return n_elems_;
    }

    /// Getter for n_points_
    inline uint n_points() const {
        return n_points_;
    }

    /// Getter for quadrature
    Quadrature *get_quadrature() const {
        return quad_;
    }

    /// Resize data tables. Method is called before reinit of patch.
    void resize_tables(uint n_elems, uint n_points) {
        n_elems_ = n_elems;
        n_points_ = n_points;
        std::vector<uint> sizes = {n_elems, n_points};
	    for (uint i=0; i<int_table_.rows(); ++i) {
	        int_table_(i) = ArenaVec<uint>(sizes[ int_sizes_[i] ], *patch_arena_);
	    }
        for (auto &elOp : operations_)
            if (elOp.size_type() != fixedSizeOp) {
                elOp.allocate_result(sizes[elOp.size_type()], *patch_arena_);
            }
        std::fill(elements_map_.begin(), elements_map_.end(), (uint)-1);
    }

    /**
     * Register element, add to coords operation
     *
     * @param coords            Coordinates of element nodes.
     * @param element_patch_idx Index of element on patch.
     */
    uint register_element(arma::mat coords, uint element_patch_idx) {
    	if (elements_map_[element_patch_idx] != (uint)-1) {
    	    // Return index of element on patch if it is registered repeatedly
    	    return elements_map_[element_patch_idx];
    	}

    	ElOp<spacedim> &op = operations_[FeBulk::BulkOps::opElCoords];
        auto coords_mat = op.result_matrix();
        std::size_t i_elem = i_elem_;
        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
                coords_mat(i_row, i_col)(i_elem) = coords(i_row, i_col);
            }

        elements_map_[element_patch_idx] = i_elem_;
        return i_elem_++;
    }

    /**
     * Register side, add to coords operations
     *
     * @param coords      Coordinates of element nodes.
     * @param side_coords Coordinates of side nodes.
     * @param side_idx    Index of side on element.
     */
    uint register_side(arma::mat elm_coords, arma::mat side_coords, uint side_idx) {
    	{
            ElOp<spacedim> &op = operations_[FeSide::SideOps::opElCoords];
            auto coords_mat = op.result_matrix();
            std::size_t i_elem = i_elem_;
            for (uint i_col=0; i_col<elm_coords.n_cols; ++i_col)
                for (uint i_row=0; i_row<elm_coords.n_rows; ++i_row) {
                    coords_mat(i_row, i_col)(i_elem) = elm_coords(i_row, i_col);
                }
    	}

    	{
            ElOp<spacedim> &op = operations_[FeSide::SideOps::opSideCoords];
            auto coords_mat = op.result_matrix();
            std::size_t i_elem = i_elem_;
            for (uint i_col=0; i_col<side_coords.n_cols; ++i_col)
                for (uint i_row=0; i_row<side_coords.n_rows; ++i_row) {
                    coords_mat(i_row, i_col)(i_elem) = side_coords(i_row, i_col);
                }
    	}

    	int_table_(3)(i_elem_) = side_idx;

        return i_elem_++;
    }

    /**
     * Register bulk point, add to int_table_
     *
     * @param elem_table_row  Index of element in temporary element table.
     * @param value_patch_idx Index of point in ElementCacheMap.
     * @param elem_idx        Index of element in Mesh.
     * @param i_point_on_elem Index of point on element
     */
    uint register_bulk_point(uint elem_table_row, uint value_patch_idx, uint elem_idx, uint i_point_on_elem) {
        uint point_pos = i_point_on_elem * n_elems_ + elem_table_row; // index of bulk point on patch
        int_table_(0)(point_pos) = value_patch_idx;
        int_table_(1)(point_pos) = elem_table_row;
        int_table_(2)(point_pos) = elem_idx;

        points_map_[value_patch_idx] = point_pos;
        return point_pos;
    }

    /**
     * Register side point, add to int_table_
     *
     * @param elem_table_row  Index of side in temporary element table.
     * @param value_patch_idx Index of point in ElementCacheMap.
     * @param elem_idx        Index of element in Mesh.
     * @param side_idx        Index of side on element.
     * @param i_point_on_side Index of point on side
     */
    uint register_side_point(uint elem_table_row, uint value_patch_idx, uint elem_idx, uint side_idx, uint i_point_on_side) {
        uint point_pos = i_point_on_side * n_elems_ + elem_table_row; // index of side point on patch
        int_table_(0)(point_pos) = value_patch_idx;
        int_table_(1)(point_pos) = elem_table_row;
        int_table_(2)(point_pos) = elem_idx;
        int_table_(4)(point_pos) = side_idx;

        points_map_[value_patch_idx] = point_pos;
        return point_pos;
    }

    /**
     * Adds accessor of new operation to operations_ vector
     *
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     * @param input_ops_vec  Indices of input operations in operations_ vector.
     * @param size_type Type of operation by size of rows
     */
    ElOp<spacedim> &make_new_op(std::initializer_list<uint> shape, ReinitFunction reinit_f, std::vector<uint> input_ops_vec, OpSizeType size_type = pointOp) {
    	ElOp<spacedim> op_accessor(this->dim_, shape, reinit_f, size_type, input_ops_vec);
    	//row_sizes_.insert(row_sizes_.end(), op_accessor.n_comp(), size_type);
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }

    /**
     * Adds accessor of new operation with fixed data size (ref data) to operations_ vector
     *
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     */
    ElOp<spacedim> &make_fixed_op(std::initializer_list<uint> shape, ReinitFunction reinit_f) {
    	return make_new_op(shape, reinit_f, {}, fixedSizeOp);
    }

    /**
     * Adds accessor of FE operation and adds operation dynamically to operations_ vector
     *
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     * @param input_ops_vec  Indices of input operations in operations_ vector.
     * @param n_dofs         Number of DOFs
     * @param size_type      Type of operation by size of rows
     */
    ElOp<spacedim> &make_fe_op(std::initializer_list<uint> shape, ReinitFunction reinit_f, std::vector<uint> input_ops_vec, uint n_dofs,
            OpSizeType size_type = pointOp) {
    	ElOp<spacedim> op_accessor(this->dim_, shape, reinit_f, size_type, input_ops_vec, n_dofs);
    	//row_sizes_.insert(row_sizes_.end(), op_accessor.n_comp() * n_dofs, size_type);
    	operations_.push_back(op_accessor);
    	return operations_[operations_.size()-1];
    }


    /**
     * Reinitializes patch data.
     *
     * Calls reinit functions defined on each operations.
     */
    void reinit_patch() {
        if (n_elems_ == 0) return; // skip if tables are empty
        for (uint i=0; i<operations_.size(); ++i)
            operations_[i].reinit_function(operations_, int_table_);
    }

    /**
     * Returns scalar output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Scalar scalar_elem_value(uint op_idx, uint point_idx) const {
        return operations_[op_idx].raw_result()(0)( int_table_(1)(points_map_[point_idx]) );
    }

    /**
     * Returns vector output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Vector vector_elem_value(uint op_idx, uint point_idx) const {
        Vector val;
        const auto &op_matrix = operations_[op_idx].raw_result();
        uint op_matrix_idx = int_table_(1)(points_map_[point_idx]);
        for (uint i=0; i<3; ++i)
            val(i) = op_matrix(i)(op_matrix_idx);
        return val;
    }

    /**
     * Returns tensor output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Tensor tensor_elem_value(uint op_idx, uint point_idx) const {
        Tensor val;
        const auto &op_matrix = operations_[op_idx].raw_result();
        uint op_matrix_idx = int_table_(1)(points_map_[point_idx]);
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = op_matrix(i+j*spacedim)(op_matrix_idx);
        return val;
    }

    /**
     * Returns scalar output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Scalar scalar_value(uint op_idx, uint point_idx, uint i_dof=0) const {
        return operations_[op_idx].raw_result()(0)(points_map_[point_idx] + i_dof*n_points_);
    }

    /**
     * Returns vector output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Vector vector_value(uint op_idx, uint point_idx, uint i_dof=0) const {
        Vector val;
        auto op_matrix = operations_[op_idx].raw_result();
        uint op_matrix_idx = points_map_[point_idx] + i_dof*n_points_;
        for (uint i=0; i<3; ++i)
            val(i) = op_matrix(i)(op_matrix_idx);
        return val;
    }

    /**
     * Returns tensor output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Tensor tensor_value(uint op_idx, uint point_idx, uint i_dof=0) const {
        Tensor val;
        auto op_matrix = operations_[op_idx].raw_result();
        uint op_matrix_idx = points_map_[point_idx] + i_dof*n_points_;
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = op_matrix(i+j*spacedim)(op_matrix_idx);
        return val;
    }

    /// Set index of FE operation in operations_ vector.
    inline void set_fe_op(FEOps fe_op, unsigned int op_vec_idx) {
        fe_ops_indices_[fe_op] = op_vec_idx;
    }

    /// Return index of FE operation in operations_ vector.
    inline unsigned int get_fe_op(FEOps fe_op) const {
        return fe_ops_indices_[fe_op];
    }

    /**
     * Performs output of data tables to stream.
     *
     * Development method.
     * @param points Allows switched off output of point table,
     * @param ints   Allows switched off output of int (connectivity to elements) table,
     */
    void print_data_tables(ostream& stream, bool points, bool ints) const {
        if (points) {
            stream << "Point vals: " << std::endl;
            for (auto &op : operations_) {
            	auto mat = op.raw_result();
                for (uint i_mat=0; i_mat<mat.rows()*mat.cols(); ++i_mat) {
                    if (mat(i_mat).data_size()==0) stream << "<empty>";
                    else {
                        const double *vals = mat(i_mat).data_ptr();
                        for (size_t i_val=0; i_val<mat(i_mat).data_size(); ++i_val)
                            stream << vals[i_val] << " ";
                    }
                    stream << std::endl;
                }
                stream << " --- end of operation ---" << std::endl;
            }
        }
        if (ints) {
            stream << "Int vals: " << int_table_.rows() << " - " << int_table_.cols() << std::endl;
            for (uint i_row=0; i_row<int_table_.rows(); ++i_row) {
                if (int_table_(i_row).data_size()==0) stream << "<empty>";
                else {
                    const uint *vals = int_table_(i_row).data_ptr();
                    for (size_t i_val=0; i_val<int_table_(i_row).data_size(); ++i_val)
                        stream << vals[i_val] << " ";
                }
                stream << std::endl;
            }
            stream << std::endl;
        }
    }

    /**
     * Performs table of fixed operations to stream.
     *
     * Development method.
     * @param bulk_side Needs set 0 (bulk) or 1 (side) for correct output of operation names.
     */
    void print_operations(ostream& stream, uint bulk_side) const {
        std::vector< std::vector<std::string> > op_names =
        {
            { "el_coords", "jacobian", "inv_jac", "jac_det", "pt_coords", "weights", "JxW", "", "", "", "", "" },
            { "el_coords", "el_jac", "el_inv_jac", "side_coords", "side_jac", "side_jac_det", "pt_coords", "weights", "JxW", "normal_vec", "", "", "", "", "" }
        };
        stream << std::setfill(' ') << " Operation" << setw(12) << "" << "Shape" << setw(2) << ""
                << "n DOFs" << setw(2) << "" << "Input operations" << endl;
        for (uint i=0; i<operations_.size(); ++i) {
            stream << " " << std::left << setw(20) << op_names[bulk_side][i] << "" << " " << setw(6) << operations_[i].format_shape() << "" << " "
                << setw(7) << operations_[i].n_dofs() << "" << " ";
            auto &input_ops = operations_[i].input_ops();
            for (auto i_o : input_ops) stream << op_names[bulk_side][i_o] << " ";
            stream << std::endl;
        }
    }

protected:
    /**
     * Hold integer values of quadrature points of defined operations.
     *
     * Table contains following rows:
     *  0: Index of quadrature point on patch
     *  1: Row of element/side in ElOp::result_ table in registration step (before expansion)
     *  2: Element idx in Mesh
     *   - last two rows are allocated only for side point table
     *  3: Index of side in element - short vector, size of column = number of sides
     *  4: Index of side in element - long vector, size of column = number of points
     * Number of used rows is given by n_points_.
     */
    IntTableArena int_table_;

    /// Set size and type of rows of int_table_, value is set implicitly in constructor of descendants
    std::vector<OpSizeType> int_sizes_;


    /// Vector of all defined operations
    std::vector<ElOp<spacedim>> operations_;

    /// Indices of FE operations in operations_ vector
    std::vector<unsigned int> fe_ops_indices_;

    uint dim_;                          ///< Dimension
    uint n_points_;                     ///< Number of points in patch
    uint n_elems_;                      ///< Number of elements in patch
    uint i_elem_;                       ///< Index of registered element in table, helper value used during patch creating.
    Quadrature *quad_;                  ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;    ///< Map of element patch indices to ElOp::result_ and int_table_ tables
    std::vector<uint> points_map_;      ///< Map of point patch indices to ElOp::result_ and int_table_ tables
    //std::vector<OpSizeType> row_sizes_; ///< hold sizes of rows by type of operation
    AssemblyArena &asm_arena_;          ///< Reference to global assembly arena of PatchFeValues
    PatchArena *patch_arena_;           ///< Pointer to global patch arena of PatchFeValues

    friend class PatchFEValues<spacedim>;
    friend class ElOp<spacedim>;
    template <class ValueType>
    friend class ElQ;
    template <class ValueType>
    friend class FeQ;
    template<unsigned int dim>
    friend class BulkValues;
    template<unsigned int dim>
    friend class SideValues;
    template<unsigned int dim>
    friend class JoinValues;
};

/**
 * @brief Class represents element or FE operations.
 */
template<unsigned int spacedim = 3>
class ElOp {
public:
    /**
     * Constructor
     *
     * Set all data members.
     */
    ElOp(uint dim, std::initializer_list<uint> shape, ReinitFunction reinit_f, OpSizeType size_type, std::vector<uint> input_ops = {}, uint n_dofs = 1)
    : dim_(dim), shape_(set_shape_vec(shape)), size_type_(size_type), input_ops_(input_ops), n_dofs_(n_dofs), reinit_func(reinit_f)
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

    /// Getter for input_ops_
    inline const std::vector<uint> &input_ops() const {
        return input_ops_;
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
        stringstream ss;
        ss << shape_[0] << "x" << shape_[1];
        return ss.str();
    }

    /// Call reinit function on element table if function is defined
    inline void reinit_function(std::vector<ElOp<spacedim>> &operations, IntTableArena &int_table) {
        reinit_func(operations, int_table);
    }

    inline void allocate_result(size_t data_size, PatchArena &arena) {
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1]);
        for (uint i=0; i<n_comp(); ++i)
            result_(i) = ArenaVec<double>(data_size, arena);
    }

    /// Return map referenced result as Eigen::Vector
    Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> result_matrix() {
        return Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>>(result_.data(), shape_[0], shape_[1]);
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
    uint dim_;                                ///< Dimension
    std::vector<uint> shape_;                 ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> result_;    ///< Result matrix of operation
    OpSizeType size_type_;                    ///< Type of operation by size of vector (element, point or fixed size)
    std::vector<uint> input_ops_;             ///< Indices of operations in PatchPointValues::operations_ vector on which ElOp is depended
    uint n_dofs_;                             ///< Number of DOFs of FE operations (or 1 in case of element operations)

    ReinitFunction reinit_func;               ///< Pointer to patch reinit function of element data table specialized by operation
};


/// Defines common functionality of reinit operations.
struct common_reinit {
	// empty base operation
	static inline void op_base(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // empty
    }

	/// Common reinit function of vector symmetric gradient on bulk and side points
    static inline void ptop_vector_sym_grad(std::vector<ElOp<3>> &operations, uint vector_sym_grad_op_idx) {
        auto &op = operations[vector_sym_grad_op_idx];
        auto grad_vector_value = operations[ op.input_ops()[0] ].result_matrix();
        auto sym_grad_value = op.result_matrix();
        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
	/// Common reinit function of vector divergence on bulk and side points
    static inline void ptop_vector_divergence(std::vector<ElOp<3>> &operations, uint vector_divergence_op_idx) {
        auto &op = operations[vector_divergence_op_idx];
        auto grad_vector_value = operations[ op.input_ops()[0] ].result_matrix();
        auto divergence_value = op.result_matrix();
        divergence_value(0,0) = grad_vector_value(0,0) + grad_vector_value(1,1) + grad_vector_value(2,2);
    }
};

/// Defines reinit operations on bulk points.
struct bulk_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_jac(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeBulk::BulkOps::opJac];
        auto jac_value = op.result_matrix();
        auto coords_value = operations[ op.input_ops()[0] ].result_matrix();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<dim; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
    template<unsigned int dim>
    static inline void elop_inv_jac(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeBulk::BulkOps::opInvJac];
        auto inv_jac_value = op.result_matrix();
        auto jac_value = operations[ op.input_ops()[0] ].result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<3, dim>(jac_value);
    }
    template<unsigned int dim>
    static inline void elop_jac_det(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeBulk::BulkOps::opJacDet];
        auto jac_det_value = op.result_matrix();
        auto jac_value = operations[ op.input_ops()[0] ].result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<3, dim>(jac_value).abs();
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_weights(std::vector<ElOp<3>> &operations, PatchArena *arena, const std::vector<double> &point_weights) {
        auto &op = operations[FeBulk::BulkOps::opWeights];
        op.allocate_result(point_weights.size(), *arena);
        auto weights_value = op.result_matrix();
        for (uint i=0; i<point_weights.size(); ++i)
            weights_value(0)(i) = point_weights[i];
    }
    static inline void ptop_JxW(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        auto &op = operations[FeBulk::BulkOps::opJxW];
        auto weights_value = operations[ op.input_ops()[0] ].result_matrix();
        auto jac_det_value = operations[ op.input_ops()[1] ].result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        auto jxw_value = op.result_matrix();
        jxw_value(0) = jxw_ovec.get_vec();
    }
    static inline void ptop_scalar_shape(std::vector<ElOp<3>> &operations,
            std::vector< std::vector<double> > shape_values, uint scalar_shape_op_idx) {
        auto &op = operations[scalar_shape_op_idx];
        uint n_dofs = shape_values.size();
        uint n_points = shape_values[0].size();
        uint n_elem = op.raw_result()(0).data_size() / n_points;

        auto shape_matrix = op.result_matrix();
        ArenaVec<double> ref_vec(n_points * n_dofs, op.raw_result()(0).arena());
        ArenaVec<double> elem_vec(n_elem, op.raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> ref_ovec(ref_vec);
        ArenaOVec<double> elem_ovec(elem_vec);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                ref_vec(i_dof * n_points + i_p) = shape_values[i_dof][i_p];
        ArenaOVec<double> shape_ovec = elem_ovec * ref_ovec;
        shape_matrix(0) = shape_ovec.get_vec();
    }
    static inline void ptop_vector_shape(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED std::vector< std::vector<arma::vec3> > shape_values, FMT_UNUSED uint vector_shape_op_idx) {
        auto &op = operations[vector_shape_op_idx];
        uint n_dofs = shape_values.size();
        uint n_points = shape_values[0].size();
        uint n_elem = op.raw_result()(0).data_size() / n_points;

        auto shape_matrix = op.result_matrix(); // result: shape 3x1
        Eigen::Matrix<ArenaVec<double>, 3, 1> ref_shape_vec;
        Eigen::Matrix<ArenaOVec<double>, 3, 1> ref_shape_ovec;
        ArenaVec<double> elem_vec(n_elem, op.raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        for (uint c=0; c<3; ++c) {
            ref_shape_vec(c) = ArenaVec<double>(n_points * n_dofs, op.raw_result()(0).arena());
            ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(c));
        }
        ArenaOVec<double> elem_ovec(elem_vec);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                for (uint c=0; c<3; ++c)
                    ref_shape_vec(c,0)(i_dof * n_points + i_p) = shape_values[i_dof][i_p](c);
        Eigen::Matrix<ArenaOVec<double>, 1, 3> shape_omatrix = elem_ovec * ref_shape_ovec.transpose();
        for (uint c=0; c<3; ++c)
            shape_matrix(c) = shape_omatrix(0,c).get_vec();
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED std::vector<ElOp<3>> &operations,
            FMT_UNUSED std::vector< std::vector<arma::vec3> > shape_values, FMT_UNUSED uint vector_shape_op_idx) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = operations[ op.input_ops()[0] ].result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    static inline void ptop_vector_piola_shape(FMT_UNUSED std::vector<ElOp<3>> &operations,
            FMT_UNUSED std::vector< std::vector<arma::vec3> > shape_values, FMT_UNUSED uint vector_shape_op_idx) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = operations[ op.input_ops()[0] ].result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(std::vector<ElOp<3>> &operations,
            std::vector< std::vector<arma::mat> > ref_shape_grads, uint scalar_shape_grads_op_idx) {
        auto &op = operations[scalar_shape_grads_op_idx];
        auto inv_jac_vec = operations[ op.input_ops()[0] ].result_matrix();
        uint n_points = ref_shape_grads.size();
        uint n_dofs = ref_shape_grads[0].size();

        Eigen::Matrix<ArenaVec<double>, dim, 1> ref_grads_vec;
        Eigen::Matrix<ArenaOVec<double>, dim, 1> ref_grads_ovec;
        for (uint i=0; i<ref_grads_vec.rows(); ++i) {
            ref_grads_vec(i) = ArenaVec<double>(n_points * n_dofs, op.raw_result()(0).arena());
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_grads_vec(i)(i_dof * n_points + i_p) = ref_shape_grads[i_p][i_dof](i);
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = op.result_matrix();
        Eigen::Matrix<ArenaOVec<double>, 3, 1> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<3; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(std::vector<ElOp<3>> &operations,
            std::vector< std::vector<arma::mat> > ref_shape_grads, uint vector_shape_grads_op_idx) {
        auto &op = operations[vector_shape_grads_op_idx];
        auto inv_jac_vec = operations[ op.input_ops()[0] ].result_matrix();
        uint n_points = ref_shape_grads.size();
        uint n_dofs = ref_shape_grads[0].size();

        Eigen::Matrix<ArenaVec<double>, dim, 3> ref_grads_vec;
        Eigen::Matrix<ArenaOVec<double>, dim, 3> ref_grads_ovec;
        for (uint i=0; i<ref_grads_vec.rows()*ref_grads_vec.cols(); ++i) {
            ref_grads_vec(i) = ArenaVec<double>(n_points * n_dofs, op.raw_result()(0).arena());
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_grads_vec(i)(i_dof * n_points + i_p) = ref_shape_grads[i_p][i_dof](i);
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = op.result_matrix();
        Eigen::Matrix<ArenaOVec<double>, 3, 3> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<9; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED std::vector<ElOp<3>> &operations,
            FMT_UNUSED std::vector< std::vector<arma::mat> > ref_shape_grads, FMT_UNUSED uint vector_shape_grads_op_idx) {
//        auto &op = operations[scalar_shape_grads_op_idx];
//        auto inv_jac_vec = operations[ op.input_ops()[0] ].result_matrix();
//        uint n_points = ref_shape_grads.size();
//        uint n_dofs = ref_shape_grads[0].size();
//
//        Eigen::Matrix<ArenaVec<double>, dim, 1> ref_grads_vec;
//        Eigen::Matrix<ArenaOVec<double>, dim, 1> ref_grads_ovec;
//        for (uint i=0; i<ref_grads_vec.rows(); ++i) {
//            ref_grads_vec(i) = ArenaVec<double>(n_points * n_dofs, op.raw_result()(0).arena());
//            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
//                for (uint i_p=0; i_p<n_points; ++i_p)
//                    ref_grads_vec(i)(i_dof * n_points + i_p) = ref_shape_grads[i_p][i_dof](i);
//            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
//        }
//
//        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
//        for (uint i=0; i<dim*3; ++i) {
//            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
//        }
//
//        auto result_vec = op.result_matrix();
//        Eigen::Matrix<ArenaOVec<double>, 3, 1> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
//        for (uint i=0; i<3; ++i) {
//            result_vec(i) = result_ovec(i).get_vec();
//        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED std::vector<ElOp<3>> &operations,
            FMT_UNUSED std::vector< std::vector<arma::mat> > ref_shape_grads, FMT_UNUSED uint vector_shape_grads_op_idx) {
//        auto &op = operations[scalar_shape_grads_op_idx];
//        auto inv_jac_vec = operations[ op.input_ops()[0] ].result_matrix();
//        uint n_points = ref_shape_grads.size();
//        uint n_dofs = ref_shape_grads[0].size();
//
//        Eigen::Matrix<ArenaVec<double>, dim, 1> ref_grads_vec;
//        Eigen::Matrix<ArenaOVec<double>, dim, 1> ref_grads_ovec;
//        for (uint i=0; i<ref_grads_vec.rows(); ++i) {
//            ref_grads_vec(i) = ArenaVec<double>(n_points * n_dofs, op.raw_result()(0).arena());
//            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
//                for (uint i_p=0; i_p<n_points; ++i_p)
//                    ref_grads_vec(i)(i_dof * n_points + i_p) = ref_shape_grads[i_p][i_dof](i);
//            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
//        }
//
//        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
//        for (uint i=0; i<dim*3; ++i) {
//            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
//        }
//
//        auto result_vec = op.result_matrix();
//        Eigen::Matrix<ArenaOVec<double>, 3, 1> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
//        for (uint i=0; i<3; ++i) {
//            result_vec(i) = result_ovec(i).get_vec();
//        }
    }
};



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_el_jac(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeSide::SideOps::opElJac];
        auto jac_value = op.result_matrix();
        auto coords_value = operations[ op.input_ops()[0] ].result_matrix();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<dim; j++)
                jac_value(i+3*j) = coords_value(i+3*(j+1)) - coords_value(i);
    }
    template<unsigned int dim>
    static inline void elop_el_inv_jac(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        auto &op = operations[FeSide::SideOps::opElInvJac];
        auto inv_jac_value = op.result_matrix();
        auto jac_value = operations[ op.input_ops()[0] ].result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<3, dim>(jac_value);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto &op = operations[FeSide::SideOps::opSideJac];
        auto jac_value = op.result_matrix();
        auto coords_value = operations[ op.input_ops()[0] ].result_matrix();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<dim-1; j++)
                jac_value(i+3*j) = coords_value(i+3*(j+1)) - coords_value(i);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac_det(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // result double, input matrix(spacedim, dim)
        auto &op = operations[FeSide::SideOps::opSideJacDet];
        auto jac_det_value = op.result_matrix();
        auto jac_value = operations[ op.input_ops()[0] ].result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<3, dim-1>(jac_value).abs();
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_weights(std::vector<ElOp<3>> &operations, PatchArena *arena, const std::vector<double> &point_weights) {
        auto &op = operations[FeSide::SideOps::opWeights];
        op.allocate_result(point_weights.size(), *arena);
        auto weights_value = op.result_matrix();
        for (uint i=0; i<point_weights.size(); ++i)
            weights_value(0)(i) = point_weights[i];
    }
    static inline void ptop_JxW(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
        auto &op = operations[FeSide::SideOps::opJxW];
        auto weights_value = operations[ op.input_ops()[0] ].result_matrix();
        auto jac_det_value = operations[ op.input_ops()[1] ].result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0,0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        auto jxw_value = op.result_matrix();
        jxw_value(0,0) = jxw_ovec.get_vec();
    }
    template<unsigned int dim>
    static inline void ptop_normal_vec(std::vector<ElOp<3>> &operations, IntTableArena &el_table) {
        auto &op = operations[FeSide::SideOps::opNormalVec];
        auto normal_value = op.result_matrix();
        auto inv_jac_value = operations[ op.input_ops()[0] ].result_matrix();
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
    static inline void ptop_scalar_shape(std::vector<ElOp<3>> &operations, IntTableArena &el_table,
            std::vector< std::vector< std::vector<double> > > shape_values, uint scalar_shape_op_idx) {
        uint n_dofs = shape_values[0][0].size();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        auto &op = operations[scalar_shape_op_idx];
        auto scalar_shape_value = op.result_matrix();
        scalar_shape_value(0) = ArenaVec<double>(n_dofs*n_patch_points, scalar_shape_value(0).arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            uint dof_shift = i_dof * n_patch_points;
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                scalar_shape_value(0)(i_pt + dof_shift) = shape_values[el_table(4)(i_pt)][i_pt / n_sides][i_dof];
        }
    }
    static inline void ptop_vector_shape(std::vector<ElOp<3>> &operations, IntTableArena &el_table,
            std::vector< std::vector< std::vector<arma::vec3> > > shape_values, uint vector_shape_op_idx) {
        uint n_dofs = shape_values[0][0].size();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        auto &op = operations[vector_shape_op_idx];
        auto vector_shape_value = op.result_matrix();
        for (uint c=0; c<3; c++)
            vector_shape_value(c) = ArenaVec<double>(n_dofs*n_patch_points, vector_shape_value(c).arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            uint dof_shift = i_dof * n_patch_points;
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                for (uint c=0; c<3; c++)
                    vector_shape_value(c)(i_pt + dof_shift) = shape_values[el_table(4)(i_pt)][i_pt / n_sides][i_dof](c);
        }
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table,
            FMT_UNUSED std::vector< std::vector< std::vector<arma::vec3> > > shape_values, FMT_UNUSED uint vector_shape_op_idx) {}
    static inline void ptop_vector_piola_shape(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table,
            FMT_UNUSED std::vector< std::vector< std::vector<arma::vec3> > > shape_values, FMT_UNUSED uint vector_shape_op_idx) {}
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(std::vector<ElOp<3>> &operations, IntTableArena &el_table,
            std::vector< std::vector< std::vector<arma::mat> > > ref_shape_grads, uint scalar_shape_grads_op_idx) {
        uint n_points = ref_shape_grads[0].size();
        uint n_dofs = ref_shape_grads[0][0].size();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Result vector
        auto &op = operations[scalar_shape_grads_op_idx];
        auto grad_scalar_shape_value = op.result_matrix();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = operations[ op.input_ops()[0] ].result_matrix();
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_dofs*n_patch_points, inv_jac_value(i).arena() );
        	for (uint j=0; j<n_dofs*n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, dim, 1> ref_shape_grads_expd;
        for (uint i=0; i<dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_dofs*n_patch_points, inv_jac_value(0).arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = (i_dof * n_points + i_pt) * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<dim; ++i_c) {
                        ref_shape_grads_expd(i_c)(i_begin + i_sd) = ref_shape_grads[el_table(3)(i_sd)][i_pt][i_dof][i_c];
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(std::vector<ElOp<3>> &operations, IntTableArena &el_table,
            std::vector< std::vector< std::vector<arma::mat> > > ref_shape_grads, uint vector_shape_grads_op_idx) {
        uint n_points = ref_shape_grads[0].size();
        uint n_dofs = ref_shape_grads[0][0].size();
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Result vector
        auto &op = operations[vector_shape_grads_op_idx];
        auto grad_scalar_shape_value = op.result_matrix();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = operations[ op.input_ops()[0] ].result_matrix();
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_dofs*n_patch_points, inv_jac_value(i).arena() );
        	for (uint j=0; j<n_dofs*n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, dim, 3> ref_shape_grads_expd;
        for (uint i=0; i<3*dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_dofs*n_patch_points, inv_jac_value(0).arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = (i_dof * n_points + i_pt) * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<3*dim; ++i_c) {
                        ref_shape_grads_expd(i_c)(i_begin + i_sd) = ref_shape_grads[el_table(3)(i_sd)][i_pt][i_dof][i_c];
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table,
            FMT_UNUSED std::vector< std::vector< std::vector<arma::mat> > > ref_shape_grads, FMT_UNUSED uint vector_shape_grads_op_idx) {}
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table,
            FMT_UNUSED std::vector< std::vector< std::vector<arma::mat> > > ref_shape_grads, FMT_UNUSED uint vector_shape_grads_op_idx) {}
};


// template specialization
template<>
inline void side_reinit::elop_sd_jac<1>(FMT_UNUSED std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
}

template<>
inline void side_reinit::elop_sd_jac_det<1>(std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
    auto &op = operations[FeSide::SideOps::opSideJacDet];
    auto result_vec = op.result_matrix();
    for (uint i=0;i<result_vec(0).data_size(); ++i) {
        result_vec(0,0)(i) = 1.0;
    }
}



namespace FeBulk {

    /// Bulk data specialization, order of item in operations_ vector corresponds to the BulkOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        /// Constructor
        PatchPointValues(uint dim, uint quad_order, AssemblyArena &asm_arena)
        : ::PatchPointValues<spacedim>(dim, asm_arena) {
            this->quad_ = new QGauss(dim, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp};
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
            /*auto &el_coords =*/ this->make_new_op( {spacedim, this->dim_+1}, &common_reinit::op_base, {}, OpSizeType::elemOp );

            /*auto &el_jac =*/ this->make_new_op( {spacedim, this->dim_}, &bulk_reinit::elop_jac<dim>, {BulkOps::opElCoords}, OpSizeType::elemOp );

            /*auto &el_inv_jac =*/ this->make_new_op( {this->dim_, spacedim}, &bulk_reinit::elop_inv_jac<dim>, {BulkOps::opJac}, OpSizeType::elemOp );

            /*auto &el_jac_det =*/ this->make_new_op( {1}, &bulk_reinit::elop_jac_det<dim>, {BulkOps::opJac}, OpSizeType::elemOp );

            // Second step: adds point values operations
            /*auto &pt_coords =*/ this->make_new_op( {spacedim}, &bulk_reinit::ptop_coords, {} );

            // use lambda reinit function
            const std::vector<double> &point_weights_vec = this->quad_->get_weights();
            auto lambda_weights = [this, point_weights_vec](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                    bulk_reinit::ptop_weights(operations, this->patch_arena_, point_weights_vec);
                };
            /*auto &weights =*/ this->make_fixed_op( {1}, lambda_weights );

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
        PatchPointValues(uint dim, uint quad_order, AssemblyArena &asm_arena)
        : ::PatchPointValues<spacedim>(dim, asm_arena) {
            this->quad_ = new QGauss(dim-1, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp, elemOp, pointOp};
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
            /*auto &el_coords =*/ this->make_new_op( {spacedim, this->dim_+1}, &common_reinit::op_base, {}, OpSizeType::elemOp );

            /*auto &el_jac =*/ this->make_new_op( {spacedim, this->dim_}, &side_reinit::elop_el_jac<dim>, {SideOps::opElCoords}, OpSizeType::elemOp );

            /*auto &el_inv_jac =*/ this->make_new_op( {this->dim_, spacedim}, &side_reinit::elop_el_inv_jac<dim>, {SideOps::opElJac}, OpSizeType::elemOp );

            // Second step: adds side values operations
            /*auto &sd_coords =*/ this->make_new_op( {spacedim, this->dim_}, &common_reinit::op_base, {}, OpSizeType::elemOp );

            /*auto &sd_jac =*/ this->make_new_op( {spacedim, this->dim_-1}, &side_reinit::elop_sd_jac<dim>, {SideOps::opSideCoords}, OpSizeType::elemOp );

            /*auto &sd_jac_det =*/ this->make_new_op( {1}, &side_reinit::elop_sd_jac_det<dim>, {SideOps::opSideJac}, OpSizeType::elemOp );

            // Third step: adds point values operations
            /*auto &coords =*/ this->make_new_op( {spacedim}, &side_reinit::ptop_coords, {} );

            // use lambda reinit function
            const std::vector<double> &point_weights_vec = this->quad_->get_weights();
            auto lambda_weights = [this, point_weights_vec](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                    side_reinit::ptop_weights(operations, this->patch_arena_, point_weights_vec);
                };
            /*auto &weights =*/ this->make_fixed_op( {1}, lambda_weights );

            /*auto &JxW =*/ this->make_new_op( {1}, &side_reinit::ptop_JxW, {SideOps::opWeights, SideOps::opSideJacDet} );

            /*auto &normal_vec =*/ this->make_new_op( {spacedim}, &side_reinit::ptop_normal_vec<dim>, {SideOps::opElInvJac}, OpSizeType::elemOp );
        }
    };

} // closing namespace FeSide



#endif /* PATCH_POINT_VALUES_HH_ */
