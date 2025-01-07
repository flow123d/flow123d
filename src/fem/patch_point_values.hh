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
template<unsigned int spacedim> class PatchOp;
template <class ValueType> class ElQ;
template <class ValueType> class FeQ;
template<unsigned int dim> class BulkValues;
template<unsigned int dim> class SideValues;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;




/// Type for conciseness
using ReinitFunction = std::function<void(PatchOp<3> *, IntTableArena &)>;


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
        /// fixed operations (reference data filled once during initialization)
        opWeights,            ///< weight of quadrature point
        opRefScalar,          ///< Scalar reference
        opRefVector,          ///< Vector reference
        opRefScalarGrad,      ///< Gradient scalar reference
        opRefVectorGrad,      ///< Gradient vector reference
        /// operations evaluated on elements
        opElCoords,           ///< coordinations of all nodes of element
        opJac,                ///< Jacobian of element
        opInvJac,             ///< inverse Jacobian
        opJacDet,             ///< determinant of Jacobian
        /// operations evaluated on quadrature points
        opCoords,             ///< coordinations of quadrature point
        opJxW,                ///< JxW value of quadrature point
        /// FE operations
        opScalarShape,        ///< Scalar shape operation
        opVectorShape,        ///< Vector shape operation
        opGradScalarShape,    ///< Scalar shape gradient
        opGradVectorShape,    ///< Vector shape gradient
        opVectorSymGrad,      ///< Vector symmetric gradient
        opVectorDivergence,   ///< Vector divergence
        opNItems              ///< Holds number of valid FE operations and value of invalid FE operation
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
        /// fixed operations (reference data filled once during initialization)
        opWeights,              ///< weight of quadrature point
        opRefScalar,            ///< Scalar reference
        opRefVector,            ///< Vector reference
        opRefScalarGrad,        ///< Gradient scalar reference
        opRefVectorGrad,        ///< Gradient vector reference
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
        opJxW,                  ///< JxW value of quadrature point
        opNormalVec,            ///< normal vector of quadrature point
        /// FE operations
        opScalarShape,         ///< Scalar shape operation
        opVectorShape,         ///< Vector shape operation
        opGradScalarShape,     ///< Scalar shape gradient
        opGradVectorShape,     ///< Vector shape gradient
        opVectorSymGrad,       ///< Vector symmetric gradient
        opVectorDivergence,    ///< Vector divergence
        opNItems               ///< Holds number of valid FE operations and value of invalid FE operation
    };
}


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
	 * Stores shared data members between PatchFeValues and PatchPoinValues
	 */
	struct PatchFeData {
	public:
	    /// Constructor
	    PatchFeData(size_t buffer_size, size_t simd_alignment)
	    : asm_arena_(buffer_size, simd_alignment), patch_arena_(nullptr) {}

	    /// Destructor
	    ~PatchFeData() {
	        if (patch_arena_!=nullptr)
	            delete patch_arena_;
	    }

	    AssemblyArena asm_arena_;    ///< Assembly arena, created and filled once during initialization
	    PatchArena *patch_arena_;    ///< Patch arena, reseted before patch reinit
	    ArenaVec<double> zero_vec_;  ///< ArenaVec of zero values of maximal length using in zero PatchPointValues construction
	};

    /**
     * Constructor
     *
     * @param dim Set dimension
     */
    PatchPointValues(uint dim, PatchFeData &patch_fe_data)
    : dim_(dim), elements_map_(300, 0), points_map_(300, 0), patch_fe_data_(patch_fe_data), needs_zero_values_(false) {
        reset();
    }

	/**
	 * Destructor.
	 */
    virtual ~PatchPointValues() {
	    for (uint i=0; i<operations_.size(); ++i)
	    	if (operations_[i] != nullptr) delete operations_[i];
	    if (needs_zero_values_)
	        delete zero_values_;
	}

    /**
     * Initialize object, set number of columns (quantities) in tables.
     */
    void initialize() {
        this->reset();
        int_table_.resize(int_sizes_.size());
        if (needs_zero_values_)
        	this->create_zero_values();
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
	        int_table_(i) = ArenaVec<uint>(sizes[ int_sizes_[i] ], *patch_fe_data_.patch_arena_);
	    }
        for (auto *elOp : operations_) {
            if (elOp == nullptr) continue;
            if (elOp->size_type() != fixedSizeOp) {
                elOp->allocate_result(sizes[elOp->size_type()], *patch_fe_data_.patch_arena_);
            }
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

    	PatchOp<spacedim> &op = *( operations_[FeBulk::BulkOps::opElCoords] );
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
            PatchOp<spacedim> &op = *( operations_[FeSide::SideOps::opElCoords] );
            auto coords_mat = op.result_matrix();
            std::size_t i_elem = i_elem_;
            for (uint i_col=0; i_col<elm_coords.n_cols; ++i_col)
                for (uint i_row=0; i_row<elm_coords.n_rows; ++i_row) {
                    coords_mat(i_row, i_col)(i_elem) = elm_coords(i_row, i_col);
                }
    	}

    	{
            PatchOp<spacedim> &op = *( operations_[FeSide::SideOps::opSideCoords] );
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
     * @param op_idx         Index of operation in operations_ vector
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     * @param size_type Type of operation by size of rows
     */
    PatchOp<spacedim> *make_new_op(uint op_idx, std::initializer_list<uint> shape, ReinitFunction reinit_f, OpSizeType size_type = pointOp) {
    	return make_fe_op(op_idx, shape, reinit_f, 1, size_type);
    }

    /**
     * Adds accessor of new operation with fixed data size (ref data) to operations_ vector
     *
     * @param op_idx         Index of operation in operations_ vector
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     */
    PatchOp<spacedim> *make_fixed_op(uint op_idx, std::initializer_list<uint> shape, ReinitFunction reinit_f) {
        return make_fe_op(op_idx, shape, reinit_f, 1, fixedSizeOp);
    }

    /**
     * Adds accessor of FE operation and adds operation dynamically to operations_ vector
     *
     * @param op_idx         Index of operation in operations_ vector
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     * @param n_dofs         Number of DOFs
     * @param size_type      Type of operation by size of rows
     */
    PatchOp<spacedim> *make_fe_op(uint op_idx, std::initializer_list<uint> shape, ReinitFunction reinit_f, uint n_dofs,
            OpSizeType size_type = pointOp) {
        if (operations_[op_idx] == nullptr) {
            std::vector<PatchOp<spacedim> *> input_ops_ptr;
            for (uint i_op : this->op_dependency_[op_idx]) {
                ASSERT_PTR(operations_[i_op]);
                input_ops_ptr.push_back(operations_[i_op]);
            }
            operations_[op_idx] = new PatchOp<spacedim>(this->dim_, shape, reinit_f, size_type, input_ops_ptr, n_dofs);
        }
    	return operations_[op_idx];
    }

    /**
     * Adds accessor of new operation with fixed data size (ref data) to operations_ vector
     *
     * @param op_idx         Index of operation in operations_ vector
     * @param shape          Shape of function output
     * @param reinit_f       Reinitialize function
     * @param n_dofs         Number of DOFs
     */
    PatchOp<spacedim> *make_fixed_fe_op(uint op_idx, std::initializer_list<uint> shape, ReinitFunction reinit_f, uint n_dofs) {
    	return make_fe_op(op_idx, shape, reinit_f, n_dofs, fixedSizeOp);
    }


    /**
     * Reinitializes patch data.
     *
     * Calls reinit functions defined on each operations.
     */
    void reinit_patch() {
        if (n_elems_ == 0) return; // skip if tables are empty
        for (uint i=0; i<operations_.size(); ++i)
            if (operations_[i] != nullptr) operations_[i]->reinit_function(operations_, int_table_);
    }

    /**
     * Returns scalar output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Scalar scalar_elem_value(uint op_idx, uint point_idx) const {
        return operations_[op_idx]->raw_result()(0)( int_table_(1)(points_map_[point_idx]) );
    }

    /**
     * Returns vector output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Vector vector_elem_value(uint op_idx, uint point_idx) const {
        Vector val;
        const auto &op_matrix = operations_[op_idx]->raw_result();
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
        const auto &op_matrix = operations_[op_idx]->raw_result();
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
        return operations_[op_idx]->raw_result()(0)(points_map_[point_idx] + i_dof*n_points_);
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
        auto op_matrix = operations_[op_idx]->raw_result();
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
        auto op_matrix = operations_[op_idx]->raw_result();
        uint op_matrix_idx = points_map_[point_idx] + i_dof*n_points_;
        for (uint i=0; i<3; ++i)
            for (uint j=0; j<3; ++j)
                val(i,j) = op_matrix(i+j*spacedim)(op_matrix_idx);
        return val;
    }

    /// return reference to assembly arena
    inline AssemblyArena &asm_arena() const {
    	return patch_fe_data_.asm_arena_;
    }

    /// return reference to patch arena
    inline PatchArena &patch_arena() const {
    	return *patch_fe_data_.patch_arena_;
    }

    /// Set flag needs_zero_values_ to true
    inline void zero_values_needed() {
        needs_zero_values_ = true;
    }


    inline PatchPointValues *zero_values() {
        ASSERT_PTR(zero_values_);
        return zero_values_;
    }

    /// Create zero_values_ object
    virtual void create_zero_values() =0;

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
            for (auto *op : operations_) {
                if (op == nullptr) continue;
            	auto mat = op->raw_result();
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
            { "weights", "ref_scalar", "ref_vector", "ref_scalar_grad", "ref_vector_grad", "el_coords", "jacobian", "inv_jac", "jac_det",
              "pt_coords", "JxW", "scalar_shape", "vector_shape", "grad_scalar_shape", "grad_vector_shape", "vector_sym_grad", "vector_divergence" },
            { "weights", "ref_scalar", "ref_vector", "ref_scalar_grad", "ref_vector_grad", "el_coords", "el_jac", "el_inv_jac", "side_coords",
              "side_jac", "side_jac_det", "pt_coords", "JxW", "normal_vec", "scalar_shape", "vector_shape", "grad_scalar_shape", "grad_vector_shape",
              "vector_sym_grad", "vector_divergence" }
        };
        stream << std::setfill(' ') << " Operation" << setw(12) << "" << "Shape" << setw(2) << ""
                << "n DOFs" << setw(2) << "" << "Input operations" << endl;
        for (uint i=0; i<operations_.size(); ++i) {
            if (operations_[i] == nullptr) continue;
            stream << " " << std::left << setw(20) << op_names[bulk_side][i] << "" << " " << setw(6) << operations_[i]->format_shape() << "" << " "
                << setw(7) << operations_[i]->n_dofs() << "" << " ";
            //auto &input_ops = operations_[i]->input_ops();
            for (auto i_o : op_dependency_[i]) stream << op_names[bulk_side][i_o] << " ";
            stream << std::endl;
        }
    }

protected:
    void create_zero_operations(std::vector<PatchOp<spacedim> *> &ref_ops);

    /**
     * Hold integer values of quadrature points of defined operations.
     *
     * Table contains following rows:
     *  0: Index of quadrature point on patch
     *  1: Row of element/side in PatchOp::result_ table in registration step (before expansion)
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
    std::vector<PatchOp<spacedim> *> operations_;

    /// Holds dependency between operations
    std::vector< std::vector<unsigned int> > op_dependency_;

    uint dim_;                          ///< Dimension
    uint n_points_;                     ///< Number of points in patch
    uint n_elems_;                      ///< Number of elements in patch
    uint i_elem_;                       ///< Index of registered element in table, helper value used during patch creating.
    Quadrature *quad_;                  ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;    ///< Map of element patch indices to PatchOp::result_ and int_table_ tables
    std::vector<uint> points_map_;      ///< Map of point patch indices to PatchOp::result_ and int_table_ tables

    PatchFeData &patch_fe_data_;        ///< Reference to PatchFeData structure shared with PatchFeValues

    bool needs_zero_values_;            ///< Flags hold whether zero_values_ object is needed
    PatchPointValues *zero_values_;     ///< PatchPointValues object returns zero values for all operations

    friend class PatchFEValues<spacedim>;
    friend class PatchOp<spacedim>;
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
class PatchOp {
public:
    /**
     * Constructor
     *
     * Set all data members.
     */
    PatchOp(uint dim, std::initializer_list<uint> shape, ReinitFunction reinit_f, OpSizeType size_type,
            std::vector<PatchOp<spacedim> *> input_ops = {}, uint n_dofs = 1)
    : dim_(dim), shape_(set_shape_vec(shape)), size_type_(size_type), input_ops_(input_ops),
      n_dofs_(n_dofs), reinit_func(reinit_f)
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

    /// Return pointer to operation of i_op index in input operation vector.
    inline PatchOp<spacedim> *input_ops(uint i_op) const {
        return input_ops_[i_op];
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
    inline void reinit_function(FMT_UNUSED std::vector<PatchOp<spacedim> *> &operations, IntTableArena &int_table) {
        reinit_func(this, int_table);
    }

    inline void allocate_result(size_t data_size, PatchArena &arena) {
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1]);
        for (uint i=0; i<n_comp(); ++i)
            result_(i) = ArenaVec<double>(data_size, arena);
    }

    inline void allocate_const_result(ArenaVec<double> &value_vec) {
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1]);
        for (uint i=0; i<n_comp(); ++i)
            result_(i) = value_vec;
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
    uint dim_;                                    ///< Dimension
    std::vector<uint> shape_;                     ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> result_;    ///< Result matrix of operation
    OpSizeType size_type_;                         ///< Type of operation by size of vector (element, point or fixed size)
    std::vector<PatchOp<spacedim> *> input_ops_;  ///< Indices of operations in PatchPointValues::operations_ vector on which PatchOp is depended
    uint n_dofs_;                                 ///< Number of DOFs of FE operations (or 1 in case of element operations)

    ReinitFunction reinit_func;                   ///< Pointer to patch reinit function of element data table specialized by operation
};


/// Defines common functionality of reinit operations.
struct common_reinit {
	// empty base operation
	static inline void op_base(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // empty
    }

    template<unsigned int dim>
    static inline void elop_jac(PatchOp<3> * result_op) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto jac_value = result_op->result_matrix();
        auto coords_value = result_op->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<dim; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }

    template<unsigned int dim>
    static inline void elop_inv_jac(PatchOp<3> * result_op) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        auto inv_jac_value = result_op->result_matrix();
        auto jac_value = result_op->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<3, dim>(jac_value);
    }

    template<unsigned int dim>
    static inline void elop_jac_det(PatchOp<3> * result_op) {
        // result double, input matrix(spacedim, dim)
        auto jac_det_value = result_op->result_matrix();
        auto jac_value = result_op->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<3, dim>(jac_value).abs();
    }

    static inline void ptop_JxW(PatchOp<3> * result_op) {
        auto weights_value = result_op->input_ops(0)->result_matrix();
        auto jac_det_value = result_op->input_ops(1)->result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0,0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        auto jxw_value = result_op->result_matrix();
        jxw_value(0,0) = jxw_ovec.get_vec();
    }

    /// Common reinit function of vector symmetric gradient on bulk and side points
    static inline void ptop_vector_sym_grad(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto grad_vector_value = result_op->input_ops(0)->result_matrix();
        auto sym_grad_value = result_op->result_matrix();
        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
	/// Common reinit function of vector divergence on bulk and side points
    static inline void ptop_vector_divergence(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto grad_vector_value = result_op->input_ops(0)->result_matrix();
        auto divergence_value = result_op->result_matrix();
        divergence_value(0,0) = grad_vector_value(0,0) + grad_vector_value(1,1) + grad_vector_value(2,2);
    }
};

/// Defines reinit operations on bulk points.
struct bulk_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_inv_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_inv_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_jac_det(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac_det<dim>(result_op);
    }

    // point operations
    static inline void ptop_coords(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_JxW(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::ptop_JxW(result_op);
    }
    static inline void ptop_scalar_shape(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto ref_vec = result_op->input_ops(0)->result_matrix()(0);

        uint n_points = ref_vec.data_size() / result_op->n_dofs(); // points per element
        uint n_elem = result_op->raw_result()(0).data_size() / n_points;

        auto shape_matrix = result_op->result_matrix();
        ArenaVec<double> elem_vec(n_elem, result_op->raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> ref_ovec(ref_vec);
        ArenaOVec<double> elem_ovec(elem_vec);
        ArenaOVec<double> shape_ovec = elem_ovec * ref_ovec;
        shape_matrix(0) = shape_ovec.get_vec();
    }
    static inline void ptop_vector_shape(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto ref_shape_vec = result_op->input_ops(0)->result_matrix();
        uint n_points = ref_shape_vec(0).data_size() / result_op->n_dofs();
        uint n_elem = result_op->raw_result()(0).data_size() / n_points;

        auto shape_matrix = result_op->result_matrix(); // result: shape 3x1
        Eigen::Matrix<ArenaOVec<double>, 3, 1> ref_shape_ovec;
        ArenaVec<double> elem_vec(n_elem, result_op->raw_result()(0).arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        for (uint c=0; c<3; ++c) {
            ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(c));
        }
        ArenaOVec<double> elem_ovec(elem_vec);
        Eigen::Matrix<ArenaOVec<double>, 1, 3> shape_omatrix = elem_ovec * ref_shape_ovec.transpose();
        for (uint c=0; c<3; ++c)
            shape_matrix(c) = shape_omatrix(0,c).get_vec();
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = op.input_ops(0).result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    static inline void ptop_vector_piola_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
//        auto &op = operations[vector_sym_grad_op_idx];
//        auto grad_vector_value = op.input_ops(0).result_matrix();
//        auto sym_grad_value = op.result_matrix();
//        sym_grad_value = 0.5 * (grad_vector_value.transpose() + grad_vector_value);
    }
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto inv_jac_vec = result_op->input_ops(0)->result_matrix();
        auto ref_grads_vec = result_op->input_ops(1)->result_matrix();

        Eigen::Matrix<ArenaOVec<double>, dim, 1> ref_grads_ovec;
        for (uint i=0; i<dim; ++i) {
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = result_op->result_matrix();
        Eigen::Matrix<ArenaOVec<double>, 3, 1> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<3; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        auto inv_jac_vec = result_op->input_ops(0)->result_matrix();
        auto ref_grads_vec = result_op->input_ops(1)->result_matrix();

        Eigen::Matrix<ArenaOVec<double>, dim, 3> ref_grads_ovec;
        for (uint i=0; i<ref_grads_vec.rows()*ref_grads_vec.cols(); ++i) {
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = result_op->result_matrix();
        Eigen::Matrix<ArenaOVec<double>, 3, 3> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<9; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
};



/// Defines reinit operations on side points.
struct side_reinit {
	// element operations
    template<unsigned int dim>
    static inline void elop_el_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_el_inv_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_inv_jac<dim>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // result matrix(spacedim, dim), input matrix(spacedim, dim+1)
        common_reinit::elop_jac<dim-1>(result_op);
    }
    template<unsigned int dim>
    static inline void elop_sd_jac_det(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::elop_jac_det<dim-1>(result_op);
    }

    // Point operations
    static inline void ptop_coords(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        // Implement
    }
    static inline void ptop_JxW(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
        common_reinit::ptop_JxW(result_op);
    }
    template<unsigned int dim>
    static inline void ptop_normal_vec(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto normal_value = result_op->result_matrix();
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();
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
    static inline void ptop_scalar_shape(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto shape_values = result_op->input_ops(0)->result_matrix();

        uint n_dofs = result_op->n_dofs();
        uint n_points = shape_values(0).data_size() / n_dofs;
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        auto scalar_shape_value = result_op->result_matrix();
        scalar_shape_value(0) = ArenaVec<double>(n_dofs*n_patch_points, scalar_shape_value(0).arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            uint dof_shift = i_dof * n_patch_points;
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                scalar_shape_value(0)(i_pt + dof_shift) = shape_values(el_table(4)(i_pt))(i_dof * n_points + i_pt / n_sides);
        }
    }
    static inline void ptop_vector_shape(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto ref_shape_vec = result_op->input_ops(0)->result_matrix();

        uint n_dofs = result_op->n_dofs();
        uint n_points = ref_shape_vec(0).data_size() / n_dofs;
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        auto vector_shape_value = result_op->result_matrix();
        for (uint c=0; c<3; c++)
            vector_shape_value(c) = ArenaVec<double>(n_dofs*n_patch_points, vector_shape_value(c).arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            uint dof_shift = i_dof * n_patch_points;
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                for (uint c=0; c<3; c++)
                    vector_shape_value(c)(i_pt + dof_shift) = ref_shape_vec(el_table(4)(i_pt),c)(i_dof * n_points + i_pt / n_sides);
        }
    }
    static inline void ptop_vector_contravariant_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    static inline void ptop_vector_piola_shape(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    template<unsigned int dim>
    static inline void ptop_scalar_shape_grads(PatchOp<3> * result_op, IntTableArena &el_table) {
        auto ref_shape_grads = result_op->input_ops(1)->result_matrix();
        auto grad_scalar_shape_value = result_op->result_matrix();  // Result vector

        uint n_dofs = result_op->n_dofs();
        uint n_points = ref_shape_grads(0).data_size() / n_dofs;
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();
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
                        ref_shape_grads_expd(i_c)(i_begin + i_sd) = ref_shape_grads(el_table(3)(i_sd),i_c)(i_dof * n_points + i_pt);
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
    template<unsigned int dim>
    static inline void ptop_vector_shape_grads(PatchOp<3> * result_op, IntTableArena &el_table) {
        // Result vector
        auto grad_scalar_shape_value = result_op->result_matrix();
        auto ref_vector_grad = result_op->input_ops(1)->result_matrix();

        uint n_dofs = result_op->n_dofs();
        uint n_points = ref_vector_grad(0).data_size() / n_dofs;
        uint n_sides = el_table(3).data_size();
        uint n_patch_points = el_table(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = result_op->input_ops(0)->result_matrix();
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
                        ref_shape_grads_expd(i_c)(i_begin + i_sd) = ref_vector_grad(el_table(3)(i_sd), i_c)(i_dof * n_points + i_pt);
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
    template<unsigned int dim>
    static inline void ptop_vector_contravariant_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
    template<unsigned int dim>
    static inline void ptop_vector_piola_shape_grads(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {}
};


// template specialization
template<>
inline void side_reinit::elop_sd_jac<1>(FMT_UNUSED PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
}

template<>
inline void side_reinit::elop_sd_jac_det<1>(PatchOp<3> * result_op, FMT_UNUSED IntTableArena &el_table) {
    auto result_vec = result_op->result_matrix();
    for (uint i=0;i<result_vec(0).data_size(); ++i) {
        result_vec(0,0)(i) = 1.0;
    }
}



namespace FeBulk {

    /**
     * Bulk data specialization, order of item in operations_ vector corresponds to the BulkOps enum
     *
     * TODO merge FeBulk::PatchPointValues and FeSide::PatchPointValues to PatchPointValues
     * when operation dependencies will be solved by DFS
     */
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        typedef typename ::PatchPointValues<spacedim>::PatchFeData PatchFeData;

        /// Constructor
        PatchPointValues(uint dim, uint quad_order, PatchFeData &patch_fe_data)
        : ::PatchPointValues<spacedim>(dim, patch_fe_data) {
            // set dependency of operations
            this->op_dependency_ = std::vector< std::vector<unsigned int> >(BulkOps::opNItems);
            this->op_dependency_[BulkOps::opJac] = {BulkOps::opElCoords};
            this->op_dependency_[BulkOps::opInvJac] = {BulkOps::opJac};
            this->op_dependency_[BulkOps::opJacDet] = {BulkOps::opJac};
            this->op_dependency_[BulkOps::opJxW] = {BulkOps::opWeights, BulkOps::opJacDet};
            this->op_dependency_[BulkOps::opScalarShape] = {BulkOps::opRefScalar};
            this->op_dependency_[BulkOps::opVectorShape] = {BulkOps::opRefVector};
            // VectorContravariant: {FeBulk::BulkOps::opRefVector, FeBulk::BulkOps::opJac}
            // VectorPiola: {FeBulk::BulkOps::opRefVector, FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet}
            this->op_dependency_[BulkOps::opGradScalarShape] = {BulkOps::opInvJac, BulkOps::opRefScalarGrad};
            this->op_dependency_[BulkOps::opGradVectorShape] = {BulkOps::opInvJac, BulkOps::opRefVectorGrad};
            // VectorContravariant: {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefVectorGrad, FeBulk::BulkOps::opJac}
            // VectorPiola: {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefVectorGrad, FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet}
            this->op_dependency_[BulkOps::opVectorSymGrad] = {BulkOps::opGradVectorShape};
            this->op_dependency_[BulkOps::opVectorDivergence] = {BulkOps::opGradVectorShape};

            this->quad_ = new QGauss(dim, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp};
            this->operations_.resize(BulkOps::opNItems, nullptr);
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

        /// Destructor
        ~PatchPointValues() {}

        /// Create zero_values_ object
        void create_zero_values() override {
		    this->zero_values_ = new PatchPointValues(this->dim_, this->operations_, this->patch_fe_data_);
        }

    private:
        /// Specialized constructor of zero values object. Do not use in other cases!
        PatchPointValues(uint dim, std::vector<PatchOp<spacedim> *> &operations, PatchFeData &patch_fe_data)
        : ::PatchPointValues<spacedim>(dim, patch_fe_data) {
            this->op_dependency_ = std::vector< std::vector<unsigned int> >(BulkOps::opNItems);
            this->create_zero_operations(operations);
        }

        /// Initialize operations vector
        template<unsigned int dim>
        void init() {
            // Fixed operation
        	auto *weights = this->make_fixed_op( BulkOps::opWeights, {1}, &common_reinit::op_base );
            // create result vector of weights operation in assembly arena
            const std::vector<double> &point_weights_vec = this->quad_->get_weights();
            weights->allocate_result(point_weights_vec.size(), this->patch_fe_data_.asm_arena_);
            auto weights_value = weights->result_matrix();
            for (uint i=0; i<point_weights_vec.size(); ++i)
                weights_value(0)(i) = point_weights_vec[i];

            // First step: adds element values operations
            /*auto &el_coords =*/ this->make_new_op( BulkOps::opElCoords, {spacedim, this->dim_+1}, &common_reinit::op_base, OpSizeType::elemOp );

            /*auto &el_jac =*/ this->make_new_op( BulkOps::opJac, {spacedim, this->dim_}, &bulk_reinit::elop_jac<dim>, OpSizeType::elemOp );

            /*auto &el_inv_jac =*/ this->make_new_op( BulkOps::opInvJac, {this->dim_, spacedim}, &bulk_reinit::elop_inv_jac<dim>, OpSizeType::elemOp );

            /*auto &el_jac_det =*/ this->make_new_op( BulkOps::opJacDet, {1}, &bulk_reinit::elop_jac_det<dim>, OpSizeType::elemOp );

            // Second step: adds point values operations
            /*auto &pt_coords =*/ this->make_new_op( BulkOps::opCoords, {spacedim}, &bulk_reinit::ptop_coords );

            /*auto &JxW =*/ this->make_new_op( BulkOps::opJxW, {1}, &bulk_reinit::ptop_JxW );
        }
    };

} // closing namespace FeBulk



namespace FeSide {

/// Bulk Side specialization, order of item in operations_ vector corresponds to the SideOps enum
    template<unsigned int spacedim = 3>
    class PatchPointValues : public ::PatchPointValues<spacedim> {
    public:
        typedef typename ::PatchPointValues<spacedim>::PatchFeData PatchFeData;

        /// Constructor
        PatchPointValues(uint dim, uint quad_order, PatchFeData &patch_fe_data)
        : ::PatchPointValues<spacedim>(dim, patch_fe_data) {
            // set dependency of operations
            this->op_dependency_ = std::vector< std::vector<unsigned int> >(SideOps::opNItems);
            this->op_dependency_[SideOps::opElJac] = {SideOps::opElCoords};
            this->op_dependency_[SideOps::opElInvJac] = {SideOps::opElJac};
            this->op_dependency_[SideOps::opSideJac] = {SideOps::opSideCoords};
            this->op_dependency_[SideOps::opSideJacDet] = {SideOps::opSideJac};
            this->op_dependency_[SideOps::opJxW] = {SideOps::opWeights, SideOps::opSideJacDet};
            this->op_dependency_[SideOps::opNormalVec] = {SideOps::opElInvJac};
            this->op_dependency_[SideOps::opScalarShape] = {SideOps::opRefScalar};
            this->op_dependency_[SideOps::opVectorShape] = {SideOps::opRefVector};
            // VectorContravariant: {FeSide::SideOps::opRefVector, FeSide::SideOps::opSideJac}
            // VectorPiola: {FeSide::SideOps::opRefVector, FeSide::SideOps::opSideJac, FeSide::SideOps::opSideJacDet}
            this->op_dependency_[SideOps::opGradScalarShape] = {SideOps::opElInvJac, SideOps::opRefScalarGrad};
            this->op_dependency_[SideOps::opGradVectorShape] = {SideOps::opElInvJac, SideOps::opRefVectorGrad};
            // VectorContravariant: {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefVectorGrad, FeSide::SideOps::opElJac}
            // VectorPiola: {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefVectorGrad, FeSide::SideOps::opElJac, FeSide::SideOps::opSideJacDet}
                // TODO ?? define and use opElJacDet
            this->op_dependency_[SideOps::opVectorSymGrad] = {SideOps::opGradVectorShape};
            this->op_dependency_[SideOps::opVectorDivergence] = {SideOps::opGradVectorShape};

            this->quad_ = new QGauss(dim-1, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp, elemOp, pointOp};
            this->operations_.resize(SideOps::opNItems, nullptr);
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

        /// Destructor
        ~PatchPointValues() {}

        /// Create zero_values_ object
        void create_zero_values() override {
        	this->zero_values_ = new PatchPointValues(this->dim_, this->operations_, this->patch_fe_data_);
        }

    private:
        /// Specialized constructor of zero values object. Do not use in other cases!
        PatchPointValues(uint dim, std::vector<PatchOp<spacedim> *> &operations, PatchFeData &patch_fe_data)
        : ::PatchPointValues<spacedim>(dim, patch_fe_data) {
            this->op_dependency_ = std::vector< std::vector<unsigned int> >(SideOps::opNItems);
            this->create_zero_operations(operations);
        }

        /// Initialize operations vector
        template<unsigned int dim>
        void init() {
            // Fixed operation
            auto *weights = this->make_fixed_op( SideOps::opWeights, {1},  &common_reinit::op_base ); //lambda_weights );
            // create result vector of weights operation in assembly arena
            const std::vector<double> &point_weights_vec = this->quad_->get_weights();
            weights->allocate_result(point_weights_vec.size(), this->patch_fe_data_.asm_arena_);
            auto weights_value = weights->result_matrix();
            for (uint i=0; i<point_weights_vec.size(); ++i)
                weights_value(0)(i) = point_weights_vec[i];

            // First step: adds element values operations
            /*auto &el_coords =*/ this->make_new_op( SideOps::opElCoords, {spacedim, this->dim_+1}, &common_reinit::op_base, OpSizeType::elemOp );

            /*auto &el_jac =*/ this->make_new_op( SideOps::opElJac, {spacedim, this->dim_}, &side_reinit::elop_el_jac<dim>, OpSizeType::elemOp );

            /*auto &el_inv_jac =*/ this->make_new_op( SideOps::opElInvJac, {this->dim_, spacedim}, &side_reinit::elop_el_inv_jac<dim>, OpSizeType::elemOp );

            // Second step: adds side values operations
            /*auto &sd_coords =*/ this->make_new_op( SideOps::opSideCoords, {spacedim, this->dim_}, &common_reinit::op_base, OpSizeType::elemOp );

            /*auto &sd_jac =*/ this->make_new_op( SideOps::opSideJac, {spacedim, this->dim_-1}, &side_reinit::elop_sd_jac<dim>, OpSizeType::elemOp );

            /*auto &sd_jac_det =*/ this->make_new_op( SideOps::opSideJacDet, {1}, &side_reinit::elop_sd_jac_det<dim>, OpSizeType::elemOp );

            // Third step: adds point values operations
            /*auto &coords =*/ this->make_new_op( SideOps::opCoords, {spacedim}, &side_reinit::ptop_coords );

            /*auto &JxW =*/ this->make_new_op( SideOps::opJxW, {1}, &side_reinit::ptop_JxW );

            /*auto &normal_vec =*/ this->make_new_op( SideOps::opNormalVec, {spacedim}, &side_reinit::ptop_normal_vec<dim>, OpSizeType::elemOp );
        }
    };

} // closing namespace FeSide



#endif /* PATCH_POINT_VALUES_HH_ */
