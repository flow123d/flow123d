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
//#include "fem/op_function.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"


template<unsigned int spacedim> class PatchOp;
template<unsigned int spacedim> class PatchFEValues;
namespace Op::Bulk::El {
    class OpCoords;
}
template <class ValueType> class ElQ;      // not necessary after merge with patch_fe_values.hh
template <class ValueType> class FeQ;
template<unsigned int dim> class BulkValues;
template<unsigned int dim> class SideValues;
using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;


/// Distinguishes operations by type and size of output rows
enum OpSizeType
{
	elemOp,      ///< operation is evaluated on elements or sides
	pointOp,     ///< operation is evaluated on quadrature points
	fixedSizeOp  ///< operation has fixed size and it is filled during initialization
};


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
    PatchPointValues(uint dim, uint quad_order, bool is_bulk, PatchFeData &patch_fe_data)
    : dim_(dim), is_bulk_(is_bulk), elements_map_(300, 0), points_map_(300, 0), patch_fe_data_(patch_fe_data),
      op_el_coords_(nullptr), op_sd_coords_(nullptr), needs_zero_values_(false) {
        reset();

        if (is_bulk) {
            switch (dim) {
            case 1:
            	init_bulk<1>(quad_order);
            	break;
            case 2:
            	init_bulk<2>(quad_order);
            	break;
            case 3:
            	init_bulk<3>(quad_order);
            	break;
            }
        } else {
            switch (dim) {
            case 1:
            	init_side<1>(quad_order);
            	break;
            case 2:
            	init_side<2>(quad_order);
            	break;
            case 3:
            	init_side<3>(quad_order);
            	break;
            }
        }
    }

	/**
	 * Destructor.
	 */
    virtual ~PatchPointValues() {
//        for (uint i=0; i<operations_.size(); ++i)
//        	if (operations_[i] != nullptr) delete operations_[i];
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
        std::vector<uint> sizes = {n_elems_, n_points_};
	    for (uint i=0; i<int_table_.rows(); ++i) {
	        int_table_(i) = ArenaVec<uint>(sizes[ int_sizes_[i] ], *patch_fe_data_.patch_arena_);
	    }
//        for (auto *elOp : operations_) {
//            if (elOp == nullptr) continue;
//            if (elOp->size_type() != fixedSizeOp) {
//                elOp->allocate_result(sizes[elOp->size_type()], *patch_fe_data_.patch_arena_);
//            }
//        }
        std::fill(elements_map_.begin(), elements_map_.end(), (uint)-1);
    }

//    /**
//     * Register element, add to coords operation
//     *
//     * @param coords            Coordinates of element nodes.
//     * @param element_patch_idx Index of element on patch.
//     */
//    uint register_element(arma::mat coords, uint element_patch_idx) {
//    	if (ppv.elements_map_[element_patch_idx] != (uint)-1) {
//    	    // Return index of element on patch if it is registered repeatedly
//    	    return ppv.elements_map_[element_patch_idx];
//    	}
//
//        auto coords_mat = ppv.op_el_coords_.result_matrix();
//        std::size_t i_elem = ppv.i_elem_;
//        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
//            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
//                coords_mat(i_row, i_col)(i_elem) = coords(i_row, i_col);
//            }
//
//        ppv.elements_map_[element_patch_idx] = ppv.i_elem_;
//        return ppv.i_elem_++;
//    }

    /**
     * Register side, add to coords operations
     *
     * @param coords      Coordinates of element nodes.
     * @param side_coords Coordinates of side nodes.
     * @param side_idx    Index of side on element.
     */
    uint register_side(FMT_UNUSED arma::mat elm_coords, FMT_UNUSED arma::mat side_coords, FMT_UNUSED uint side_idx) {
//    	{
//            PatchOp<spacedim> &op = *( operations_[FeSide::SideOps::opElCoords] );
//            auto coords_mat = op.result_matrix();
//            std::size_t i_elem = i_elem_;
//            for (uint i_col=0; i_col<elm_coords.n_cols; ++i_col)
//                for (uint i_row=0; i_row<elm_coords.n_rows; ++i_row) {
//                    coords_mat(i_row, i_col)(i_elem) = elm_coords(i_row, i_col);
//                }
//    	}
//
//    	{
//            PatchOp<spacedim> &op = *( operations_[FeSide::SideOps::opSideCoords] );
//            auto coords_mat = op.result_matrix();
//            std::size_t i_elem = i_elem_;
//            for (uint i_col=0; i_col<side_coords.n_cols; ++i_col)
//                for (uint i_row=0; i_row<side_coords.n_rows; ++i_row) {
//                    coords_mat(i_row, i_col)(i_elem) = side_coords(i_row, i_col);
//                }
//    	}
//
//    	int_table_(3)(i_elem_) = side_idx;

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
    PatchOp<spacedim> *make_fe_op(FMT_UNUSED uint op_idx, FMT_UNUSED std::initializer_list<uint> shape, FMT_UNUSED ReinitFunction reinit_f, FMT_UNUSED uint n_dofs,
            FMT_UNUSED OpSizeType size_type = pointOp) {
//        if (operations_[op_idx] == nullptr) {
//            std::vector<PatchOp<spacedim> *> input_ops_ptr;
//            for (uint i_op : this->op_dependency_[op_idx]) {
//                ASSERT_PTR(operations_[i_op]);
//                input_ops_ptr.push_back(operations_[i_op]);
//            }
//            operations_[op_idx] = new PatchOp<spacedim>(this->dim_, shape, reinit_f, size_type, input_ops_ptr, n_dofs);
//        }
//    	return operations_[op_idx];
    	return nullptr;
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
//        if (n_elems_ == 0) return; // skip if tables are empty
//        for (uint i=0; i<operations_.size(); ++i) {
//            if (operations_[i] != nullptr) operations_[i]->reinit_function(operations_, int_table_);
//        }
    }

    /**
     * Returns scalar output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Scalar scalar_elem_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx) const {
//        return operations_[op_idx]->raw_result()(0)( int_table_(1)(points_map_[point_idx]) );
        return 0.0;
    }

    /**
     * Returns vector output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Vector vector_elem_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx) const {
        Vector val;
//        const auto &op_matrix = operations_[op_idx]->raw_result();
//        uint op_matrix_idx = int_table_(1)(points_map_[point_idx]);
//        for (uint i=0; i<3; ++i)
//            val(i) = op_matrix(i)(op_matrix_idx);
        return val;
    }

    /**
     * Returns tensor output value of data stored by elements.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    inline Tensor tensor_elem_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx) const {
        Tensor val;
//        const auto &op_matrix = operations_[op_idx]->raw_result();
//        uint op_matrix_idx = int_table_(1)(points_map_[point_idx]);
//        for (uint i=0; i<3; ++i)
//            for (uint j=0; j<3; ++j)
//                val(i,j) = op_matrix(i+j*spacedim)(op_matrix_idx);
        return val;
    }

    /**
     * Returns scalar output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Scalar scalar_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx, FMT_UNUSED uint i_dof=0) const {
//        return operations_[op_idx]->raw_result()(i_dof)(points_map_[point_idx]);
        return 0.0;
    }

    /**
     * Returns vector output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Vector vector_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx, FMT_UNUSED uint i_dof=0) const {
        Vector val;
//        auto op_matrix = operations_[op_idx]->raw_result();
//        uint op_matrix_idx = points_map_[point_idx];
//        for (uint i=0; i<3; ++i)
//            val(i) = op_matrix(i + 3*i_dof)(op_matrix_idx);
        return val;
    }

    /**
     * Returns tensor output value on point.
     *
     * @param op_idx      Index of operation in operations vector
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    inline Tensor tensor_value(FMT_UNUSED uint op_idx, FMT_UNUSED uint point_idx, FMT_UNUSED uint i_dof=0) const {
        Tensor val;
//        auto op_matrix = operations_[op_idx]->raw_result();
//        uint op_matrix_idx = points_map_[point_idx];
//        for (uint i=0; i<9; ++i)
//            val(i) = op_matrix(i+9*i_dof)(op_matrix_idx);
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
    void create_zero_values() {
	    zero_values_ = new PatchPointValues(dim_, operations_, is_bulk_, patch_fe_data_);
    }

    /**
     * Performs output of data tables to stream.
     *
     * Development method.
     * @param points Allows switched off output of point table,
     * @param ints   Allows switched off output of int (connectivity to elements) table,
     */
    void print_data_tables(FMT_UNUSED std::ostream& stream, FMT_UNUSED bool points, FMT_UNUSED bool ints) const {
//        if (points) {
//            stream << "Point vals: " << std::endl;
//            for (auto *op : operations_) {
//                if (op == nullptr) continue;
//            	auto mat = op->raw_result();
//                for (uint i_mat=0; i_mat<mat.rows()*mat.cols(); ++i_mat) {
//                    if (mat(i_mat).data_size()==0) stream << "<empty>";
//                    else {
//                        const double *vals = mat(i_mat).data_ptr();
//                        for (size_t i_val=0; i_val<mat(i_mat).data_size(); ++i_val)
//                            stream << vals[i_val] << " ";
//                    }
//                    stream << std::endl;
//                }
//                stream << " --- end of operation ---" << std::endl;
//            }
//        }
//        if (ints) {
//            stream << "Int vals: " << int_table_.rows() << " - " << int_table_.cols() << std::endl;
//            for (uint i_row=0; i_row<int_table_.rows(); ++i_row) {
//                if (int_table_(i_row).data_size()==0) stream << "<empty>";
//                else {
//                    const uint *vals = int_table_(i_row).data_ptr();
//                    for (size_t i_val=0; i_val<int_table_(i_row).data_size(); ++i_val)
//                        stream << vals[i_val] << " ";
//                }
//                stream << std::endl;
//            }
//            stream << std::endl;
//        }
    }

    /**
     * Performs table of fixed operations to stream.
     *
     * Development method.
     * @param bulk_side Needs set 0 (bulk) or 1 (side) for correct output of operation names.
     */
    void print_operations(FMT_UNUSED std::ostream& stream, FMT_UNUSED uint bulk_side) const {
//        std::vector< std::vector<std::string> > op_names =
//        {
//            { "weights", "ref_scalar", "ref_vector", "ref_scalar_grad", "ref_vector_grad", "el_coords", "jacobian", "inv_jac", "jac_det",
//              "pt_coords", "JxW", "scalar_shape", "vector_shape", "grad_scalar_shape", "grad_vector_shape", "vector_sym_grad", "vector_divergence" },
//            { "weights", "ref_scalar", "ref_vector", "ref_scalar_grad", "ref_vector_grad", "el_coords", "el_jac", "el_inv_jac", "side_coords",
//              "side_jac", "side_jac_det", "pt_coords", "JxW", "normal_vec", "scalar_shape", "vector_shape", "grad_scalar_shape", "grad_vector_shape",
//              "vector_sym_grad", "vector_divergence" }
//        };
//        stream << std::setfill(' ') << " Operation" << std::setw(12) << "" << "Shape" << std::setw(2) << ""
//                << "n DOFs" << std::setw(2) << "" << "Input operations" << std::endl;
//        for (uint i=0; i<operations_.size(); ++i) {
//            if (operations_[i] == nullptr) continue;
//            stream << " " << std::left << std::setw(20) << op_names[bulk_side][i] << "" << " " << std::setw(6) << operations_[i]->format_shape() << "" << " "
//                << std::setw(7) << operations_[i]->n_dofs() << "" << " ";
//            //auto &input_ops = operations_[i]->input_ops();
//            for (auto i_o : op_dependency_[i]) stream << op_names[bulk_side][i_o] << " ";
//            stream << std::endl;
//        }
    }

protected:
    /// Specialized constructor of zero values object. Do not use in other cases!
    PatchPointValues(uint dim, FMT_UNUSED std::vector<PatchOp<spacedim> *> &ref_ops, bool is_bulk, PatchFeData &patch_fe_data)
    : dim_(dim), is_bulk_(is_bulk), elements_map_(300, 0), points_map_(300, 0), patch_fe_data_(patch_fe_data),
      op_el_coords_(nullptr), op_sd_coords_(nullptr), needs_zero_values_(false) {
        reset();

        if (is_bulk) op_dependency_ = std::vector< std::vector<unsigned int> >(FeBulk::BulkOps::opNItems);
        else op_dependency_ = std::vector< std::vector<unsigned int> >(FeSide::SideOps::opNItems);

//        operations_.resize(ref_ops.size(), nullptr);
//        for (uint i_op = 0; i_op < ref_ops.size(); ++i_op ) {
//            auto *op = ref_ops[i_op];
//            if (op == nullptr) continue;
//
//            auto *new_op = make_fe_op(i_op, {op->shape()[0], op->shape()[1]}, &common_reinit::op_base, op->n_dofs(), op->size_type());
//            new_op->allocate_const_result(patch_fe_data_.zero_vec_);
//        }
    }

    /// Specialized initialization part of bulk PatchPointValues
    template<unsigned int dim>
    void init_bulk(uint quad_order) {
        // set dependency of operations
        this->op_dependency_ = std::vector< std::vector<unsigned int> >(FeBulk::BulkOps::opNItems);
        this->op_dependency_[FeBulk::BulkOps::opJac] = {FeBulk::BulkOps::opElCoords};
        this->op_dependency_[FeBulk::BulkOps::opInvJac] = {FeBulk::BulkOps::opJac};
        this->op_dependency_[FeBulk::BulkOps::opJacDet] = {FeBulk::BulkOps::opJac};
        this->op_dependency_[FeBulk::BulkOps::opJxW] = {FeBulk::BulkOps::opWeights, FeBulk::BulkOps::opJacDet};
        this->op_dependency_[FeBulk::BulkOps::opScalarShape] = {FeBulk::BulkOps::opRefScalar};
        this->op_dependency_[FeBulk::BulkOps::opVectorShape] = {FeBulk::BulkOps::opRefVector};
        // VectorContravariant: {FeBulk::BulkOps::opRefVector, FeBulk::BulkOps::opJac}
        // VectorPiola: {FeBulk::BulkOps::opRefVector, FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet}
        this->op_dependency_[FeBulk::BulkOps::opGradScalarShape] = {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefScalarGrad};
        this->op_dependency_[FeBulk::BulkOps::opGradVectorShape] = {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefVectorGrad};
        // VectorContravariant: {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefVectorGrad, FeBulk::BulkOps::opJac}
        // VectorPiola: {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opRefVectorGrad, FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet}
        this->op_dependency_[FeBulk::BulkOps::opVectorSymGrad] = {FeBulk::BulkOps::opGradVectorShape};
        this->op_dependency_[FeBulk::BulkOps::opVectorDivergence] = {FeBulk::BulkOps::opGradVectorShape};

        // initialize vectors
        this->quad_ = new QGauss(dim, 2*quad_order);
        this->int_sizes_ = {pointOp, pointOp, pointOp};
//        this->operations_.resize(FeBulk::BulkOps::opNItems, nullptr);
//
//        // Create fixed operation
//    	auto *weights = this->make_fixed_op( FeBulk::BulkOps::opWeights, {1}, &common_reinit::op_base );
//        // create result vector of weights operation in assembly arena
//        const std::vector<double> &point_weights_vec = this->quad_->get_weights();
//        weights->allocate_result(point_weights_vec.size(), this->patch_fe_data_.asm_arena_);
//        auto weights_value = weights->result_matrix();
//        for (uint i=0; i<point_weights_vec.size(); ++i)
//            weights_value(0)(i) = point_weights_vec[i];
//
//        // First step: adds element values operations
//        /*auto &el_coords =*/ this->make_new_op( FeBulk::BulkOps::opElCoords, {spacedim, this->dim_+1}, &common_reinit::op_base, OpSizeType::elemOp );
//
//        /*auto &el_jac =*/ this->make_new_op( FeBulk::BulkOps::opJac, {spacedim, this->dim_}, &bulk_reinit::elop_jac<dim>, OpSizeType::elemOp );
//
//        /*auto &el_inv_jac =*/ this->make_new_op( FeBulk::BulkOps::opInvJac, {this->dim_, spacedim}, &bulk_reinit::elop_inv_jac<dim>, OpSizeType::elemOp );
//
//        /*auto &el_jac_det =*/ this->make_new_op( FeBulk::BulkOps::opJacDet, {1}, &bulk_reinit::elop_jac_det<dim>, OpSizeType::elemOp );
//
//        // Second step: adds point values operations
//        /*auto &pt_coords =*/ this->make_new_op( FeBulk::BulkOps::opCoords, {spacedim}, &bulk_reinit::ptop_coords );
//
//        /*auto &JxW =*/ this->make_new_op( FeBulk::BulkOps::opJxW, {1}, &bulk_reinit::ptop_JxW );
    }

    /// Specialized initialization part of side PatchPointValues
    template<unsigned int dim>
    void init_side(uint quad_order) {
        // set dependency of operations
        this->op_dependency_ = std::vector< std::vector<unsigned int> >(FeSide::SideOps::opNItems);
        this->op_dependency_[FeSide::SideOps::opElJac] = {FeSide::SideOps::opElCoords};
        this->op_dependency_[FeSide::SideOps::opElInvJac] = {FeSide::SideOps::opElJac};
        this->op_dependency_[FeSide::SideOps::opSideJac] = {FeSide::SideOps::opSideCoords};
        this->op_dependency_[FeSide::SideOps::opSideJacDet] = {FeSide::SideOps::opSideJac};
        this->op_dependency_[FeSide::SideOps::opJxW] = {FeSide::SideOps::opWeights, FeSide::SideOps::opSideJacDet};
        this->op_dependency_[FeSide::SideOps::opNormalVec] = {FeSide::SideOps::opElInvJac};
        this->op_dependency_[FeSide::SideOps::opScalarShape] = {FeSide::SideOps::opRefScalar};
        this->op_dependency_[FeSide::SideOps::opVectorShape] = {FeSide::SideOps::opRefVector};
        // VectorContravariant: {FeSide::SideOps::opRefVector, FeSide::SideOps::opSideJac}
        // VectorPiola: {FeSide::SideOps::opRefVector, FeSide::SideOps::opSideJac, FeSide::SideOps::opSideJacDet}
        this->op_dependency_[FeSide::SideOps::opGradScalarShape] = {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefScalarGrad};
        this->op_dependency_[FeSide::SideOps::opGradVectorShape] = {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefVectorGrad};
        // VectorContravariant: {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefVectorGrad, FeSide::SideOps::opElJac}
        // VectorPiola: {FeSide::SideOps::opElInvJac, FeSide::SideOps::opRefVectorGrad, FeSide::SideOps::opElJac, FeSide::SideOps::opSideJacDet}
            // TODO ?? define and use opElJacDet
        this->op_dependency_[FeSide::SideOps::opVectorSymGrad] = {FeSide::SideOps::opGradVectorShape};
        this->op_dependency_[FeSide::SideOps::opVectorDivergence] = {FeSide::SideOps::opGradVectorShape};

        // initialize vectors
        this->quad_ = new QGauss(dim-1, 2*quad_order);
    	this->int_sizes_ = {pointOp, pointOp, pointOp, elemOp, pointOp};
//        this->operations_.resize(FeSide::SideOps::opNItems, nullptr);
//
//        // Create fixed operation
//        auto *weights = this->make_fixed_op( FeSide::SideOps::opWeights, {1},  &common_reinit::op_base ); //lambda_weights );
//        // create result vector of weights operation in assembly arena
//        const std::vector<double> &point_weights_vec = this->quad_->get_weights();
//        weights->allocate_result(point_weights_vec.size(), this->patch_fe_data_.asm_arena_);
//        auto weights_value = weights->result_matrix();
//        for (uint i=0; i<point_weights_vec.size(); ++i)
//            weights_value(0)(i) = point_weights_vec[i];
//
//        // First step: adds element values operations
//        /*auto &el_coords =*/ this->make_new_op( FeSide::SideOps::opElCoords, {spacedim, this->dim_+1}, &common_reinit::op_base, OpSizeType::elemOp );
//
//        /*auto &el_jac =*/ this->make_new_op( FeSide::SideOps::opElJac, {spacedim, this->dim_}, &side_reinit::elop_el_jac<dim>, OpSizeType::elemOp );
//
//        /*auto &el_inv_jac =*/ this->make_new_op( FeSide::SideOps::opElInvJac, {this->dim_, spacedim}, &side_reinit::elop_el_inv_jac<dim>, OpSizeType::elemOp );
//
//        // Second step: adds side values operations
//        /*auto &sd_coords =*/ this->make_new_op( FeSide::SideOps::opSideCoords, {spacedim, this->dim_}, &common_reinit::op_base, OpSizeType::elemOp );
//
//        /*auto &sd_jac =*/ this->make_new_op( FeSide::SideOps::opSideJac, {spacedim, this->dim_-1}, &side_reinit::elop_sd_jac<dim>, OpSizeType::elemOp );
//
//        /*auto &sd_jac_det =*/ this->make_new_op( FeSide::SideOps::opSideJacDet, {1}, &side_reinit::elop_sd_jac_det<dim>, OpSizeType::elemOp );
//
//        // Third step: adds point values operations
//        /*auto &coords =*/ this->make_new_op( FeSide::SideOps::opCoords, {spacedim}, &side_reinit::ptop_coords );
//
//        /*auto &JxW =*/ this->make_new_op( FeSide::SideOps::opJxW, {1}, &side_reinit::ptop_JxW );
//
//        /*auto &normal_vec =*/ this->make_new_op( FeSide::SideOps::opNormalVec, {spacedim}, &side_reinit::ptop_normal_vec<dim>, OpSizeType::elemOp );
    }

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
    bool is_bulk_;                      ///< Flag, bulk or side PatchPointValues
    uint n_points_;                     ///< Number of points in patch
    uint n_elems_;                      ///< Number of elements in patch
    uint i_elem_;                       ///< Index of registered element in table, helper value used during patch creating.
    Quadrature *quad_;                  ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;    ///< Map of element patch indices to PatchOp::result_ and int_table_ tables
    std::vector<uint> points_map_;      ///< Map of point patch indices to PatchOp::result_ and int_table_ tables

    PatchFeData &patch_fe_data_;        ///< Reference to PatchFeData structure shared with PatchFeValues

    PatchOp<spacedim> *op_el_coords_;   ///< Pointer to element coords operations (used during patch reinit)
	PatchOp<spacedim> *op_sd_coords_;   ///< Pointer to side coords operations (used during patch reinit)

    bool needs_zero_values_;            ///< Flags hold whether zero_values_ object is needed
    PatchPointValues *zero_values_;     ///< PatchPointValues object returns zero values for all operations

    friend class PatchFEValues<spacedim>;
    friend class PatchOp<spacedim>;
    friend class Op::Bulk::El::OpCoords;
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


#endif /* PATCH_POINT_VALUES_HH_ */
