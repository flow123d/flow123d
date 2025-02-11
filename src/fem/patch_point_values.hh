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
#include "quadrature/quadrature_lib.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"


template<unsigned int spacedim> class PatchOp;
template<unsigned int spacedim> class PatchFEValues;
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
      needs_zero_values_(false) {
        reset();

        if (is_bulk) {
            this->quad_ = new QGauss(dim, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp};
        } else {
            this->quad_ = new QGauss(dim-1, 2*quad_order);
            this->int_sizes_ = {pointOp, pointOp, pointOp, elemOp, pointOp};
        }
    }

	/**
	 * Destructor.
	 */
    virtual ~PatchPointValues() {
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
        elem_list_.clear();
        side_list_.clear();
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
        std::fill(elements_map_.begin(), elements_map_.end(), (uint)-1);
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
	    ASSERT_PERMANENT(false);
        // zero_values_ = new PatchPointValues(dim_, operations_, is_bulk_, patch_fe_data_);
    }

//protected:
    /// Specialized constructor of zero values object. Do not use in other cases!
    PatchPointValues(uint dim, FMT_UNUSED std::vector<PatchOp<spacedim> *> &ref_ops, bool is_bulk, PatchFeData &patch_fe_data)
    : dim_(dim), is_bulk_(is_bulk), elements_map_(300, 0), points_map_(300, 0), patch_fe_data_(patch_fe_data),
      needs_zero_values_(false) {
        reset();

//        if (is_bulk) op_dependency_ = std::vector< std::vector<unsigned int> >(FeBulk::BulkOps::opNItems);
//        else op_dependency_ = std::vector< std::vector<unsigned int> >(FeSide::SideOps::opNItems);
//
//        operations_.resize(ref_ops.size(), nullptr);
//        for (uint i_op = 0; i_op < ref_ops.size(); ++i_op ) {
//            auto *op = ref_ops[i_op];
//            if (op == nullptr) continue;
//
//            auto *new_op = make_fe_op(i_op, {op->shape()[0], op->shape()[1]}, &common_reinit::op_base, op->n_dofs(), op->size_type());
//            new_op->allocate_const_result(patch_fe_data_.zero_vec_);
//        }
    }

    /// Dependencies of bulk operations
/*      this->op_dependency_ = std::vector< std::vector<unsigned int> >(FeBulk::BulkOps::opNItems);
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
*/

    /// Dependencies of side operations -
/*      this->op_dependency_ = std::vector< std::vector<unsigned int> >(FeSide::SideOps::opNItems);
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
        this->operations_.resize(FeSide::SideOps::opNItems, nullptr);

        // Create fixed operation
        auto *weights = this->make_fixed_op( FeSide::SideOps::opWeights, {1},  &common_reinit::op_base ); //lambda_weights );
        // create result vector of weights operation in assembly arena
        const std::vector<double> &point_weights_vec = this->quad_->get_weights();
        weights->allocate_result(point_weights_vec.size(), this->patch_fe_data_.asm_arena_);
        auto weights_value = weights->result_matrix();
        for (uint i=0; i<point_weights_vec.size(); ++i)
            weights_value(0)(i) = point_weights_vec[i];

        // First step: adds element values operations
        auto &el_coords = this->make_new_op( FeSide::SideOps::opElCoords, {spacedim, this->dim_+1}, &common_reinit::op_base, OpSizeType::elemOp );

        auto &el_jac = this->make_new_op( FeSide::SideOps::opElJac, {spacedim, this->dim_}, &side_reinit::elop_el_jac<dim>, OpSizeType::elemOp );

        auto &el_inv_jac / this->make_new_op( FeSide::SideOps::opElInvJac, {this->dim_, spacedim}, &side_reinit::elop_el_inv_jac<dim>, OpSizeType::elemOp );

        // Second step: adds side values operations
        auto &sd_coords = this->make_new_op( FeSide::SideOps::opSideCoords, {spacedim, this->dim_}, &common_reinit::op_base, OpSizeType::elemOp );

        auto &sd_jac = this->make_new_op( FeSide::SideOps::opSideJac, {spacedim, this->dim_-1}, &side_reinit::elop_sd_jac<dim>, OpSizeType::elemOp );

        auto &sd_jac_det = this->make_new_op( FeSide::SideOps::opSideJacDet, {1}, &side_reinit::elop_sd_jac_det<dim>, OpSizeType::elemOp );

        // Third step: adds point values operations
        auto &coords = this->make_new_op( FeSide::SideOps::opCoords, {spacedim}, &side_reinit::ptop_coords );

        auto &JxW = this->make_new_op( FeSide::SideOps::opJxW, {1}, &side_reinit::ptop_JxW );

        auto &normal_vec = this->make_new_op( FeSide::SideOps::opNormalVec, {spacedim}, &side_reinit::ptop_normal_vec<dim>, OpSizeType::elemOp );
*/


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
    //std::vector<PatchOp<spacedim> *> operations_;

    uint dim_;                          ///< Dimension
    bool is_bulk_;                      ///< Flag, bulk or side PatchPointValues
    uint n_points_;                     ///< Number of points in patch
    uint n_elems_;                      ///< Number of elements in patch
    uint i_elem_;                       ///< Index of registered element in table, helper value used during patch creating.
    Quadrature *quad_;                  ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;    ///< Map of element patch indices to PatchOp::result_ and int_table_ tables
    std::vector<uint> points_map_;      ///< Map of point patch indices to PatchOp::result_ and int_table_ tables

    PatchFeData &patch_fe_data_;        ///< Reference to PatchFeData structure shared with PatchFeValues

	std::vector<ElementAccessor<3>> elem_list_; ///< List of elements on patch
	std::vector<Side> side_list_;               ///< List of sides on patch

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


#endif /* PATCH_POINT_VALUES_HH_ */
