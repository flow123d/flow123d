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
    : elements_map_(300, 0), points_map_(300, 0), patch_fe_data_(patch_fe_data) {
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
    virtual ~PatchPointValues() {}

    /**
     * Initialize object, set number of columns (quantities) in tables.
     */
    void initialize() {
        this->reset();
        int_table_.resize(int_sizes_.size());
    }

    /// Reset number of columns (points and elements)
    inline void reset() {
        n_points_ = 0;
        n_elems_ = 0;
        i_elem_ = 0;
        elem_list_.clear();
        side_list_.clear();
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

    /// return reference to assembly arena
    inline AssemblyArena &asm_arena() const {
    	return patch_fe_data_.asm_arena_;
    }

    /// return reference to patch arena
    inline PatchArena &patch_arena() const {
    	return *patch_fe_data_.patch_arena_;
    }

//protected:

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


    uint n_points_;                     ///< Number of points in patch
    uint n_elems_;                      ///< Number of elements in patch
    uint i_elem_;                       ///< Index of registered element in table, helper value used during patch creating.
    Quadrature *quad_;                  ///< Quadrature of given dimension and order passed in constructor.

    std::vector<uint> elements_map_;    ///< Map of element patch indices to PatchOp::result_ and int_table_ tables
    std::vector<uint> points_map_;      ///< Map of point patch indices to PatchOp::result_ and int_table_ tables

    PatchFeData &patch_fe_data_;        ///< Reference to PatchFeData structure shared with PatchFeValues

	std::vector<ElementAccessor<3>> elem_list_; ///< List of elements on patch
	std::vector<Side> side_list_;               ///< List of sides on patch

    friend class PatchFEValues<spacedim>;
    friend class PatchOp<spacedim>;
};


#endif /* PATCH_POINT_VALUES_HH_ */
