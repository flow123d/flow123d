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
#include <unordered_set>

#include "fem/eigen_tools.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"
#include "mesh/accessors.hh"



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



struct RevertibleValue {
public:
    /// Default constructor
    RevertibleValue() : permanent_(0), temporary_(0) {}

    /// Copy constructor
    RevertibleValue(const RevertibleValue &other)
    : permanent_(other.permanent_), temporary_(other.temporary_) {}


    /// Declaration of operators
    inline RevertibleValue &operator= (RevertibleValue &other) {
    	permanent_ = other.permanent_;
    	temporary_ = other.temporary_;
        return *this;
    }

    inline RevertibleValue& operator++ ()
    {
    	temporary_++;
        return (*this);
    }

    inline RevertibleValue& operator+= (std::size_t inc_val)
    {
    	temporary_ += inc_val;
        return (*this);
    }

    inline std::size_t operator() () const
    {
        return permanent_;
    }


    /// Reset value to zero
    void reset() {
        permanent_ = 0;
        temporary_ = 0;
    }

    /// Revert temporary value.
    inline void revert_temporary() {
        temporary_ = permanent_;
    }

    /// Finalize temporary value.
    inline void make_permanent() {
        permanent_ = temporary_;
    }

    /// Return temporary value.
    inline std::size_t temporary_value() const
    {
        return temporary_;
    }

private:
    std::size_t permanent_;
    std::size_t temporary_;
};


/**
 * Data class, holds common data of elements of one dimension of bulk and side PatchPointValues.
 */
template<unsigned int spacedim = 3>
class ElemsDimOata {
public:
    ElemsDimOata()
    : i_elem_(0) {}

    /// Reset data
    void reset() {
        n_elems_.clear();
        i_elem_ = 0;
        elem_list_.clear();
    }

    std::unordered_set<uint> n_elems_;                  ///< Holds number of elements of one dim n patch, unordered_set ensures control of duplicity.
    uint i_elem_;                                       ///< Index of registered element in table, helper value used during patch creating
    std::vector<ElementAccessor<spacedim>> elem_list_;  ///< List of elements on patch
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
    PatchPointValues(ElemsDimOata<spacedim> *elems_dim_data, bool is_bulk)
    : elems_dim_data_(elems_dim_data), points_map_(300, 0) {
        reset();

        if (is_bulk) {
            this->int_sizes_ = {pointOp, pointOp, pointOp};
        } else {
            this->int_sizes_ = {pointOp, pointOp, pointOp, elemOp, pointOp, elemOp};
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
        n_points_.reset();
        n_sides_.reset();
        i_side_ = 0;
        elems_dim_data_->reset();
        side_list_.clear();
    }

    /// Getter for n_elems_
    inline uint n_elems() const {
        return elems_dim_data_->n_elems_.size();
    }

    /// Getter for n_sides_
    inline uint n_sides() const {
        return n_sides_();
    }

    /// Getter for n_points_
    inline uint n_points() const {
        return n_points_();
    }

    /// Resize data tables. Method is called before reinit of patch.
    void resize_tables(PatchArena &patch_arena) {
        std::vector<std::size_t> sizes = {n_sides_(), n_points_()};
	    for (uint i=0; i<int_table_.rows(); ++i) {
	        int_table_(i) = ArenaVec<uint>(sizes[ int_sizes_[i] ], patch_arena);
	    }
    }

    /**
     * Register bulk point, add to int_table_
     *
     * @param patch_elm_idx     Index of element on patch (in int_table_).
     * @param elm_cache_map_idx Index of point in ElementCacheMap.
     * @param elem_idx          Index of element in Mesh.
     * @param i_point_on_elem   Index of point on element
     */
    uint register_bulk_point(uint patch_elm_idx, uint elm_cache_map_idx, uint elem_idx, uint i_point_on_elem) {
        uint point_pos = i_point_on_elem * n_elems() + patch_elm_idx; // index of bulk point on patch
        int_table_(0)(point_pos) = elm_cache_map_idx;
        int_table_(1)(point_pos) = patch_elm_idx;
        int_table_(2)(point_pos) = elem_idx;

        points_map_[elm_cache_map_idx] = point_pos;
        return point_pos;
    }

    /**
     * Register side point, add to int_table_
     *
     * @param patch_side_idx     Index of side on patch (in int_table_).
     * @param elm_cache_map_idx  Index of point in ElementCacheMap.
     * @param elem_idx           Index of element in Mesh.
     * @param side_idx           Index of side on element.
     * @param i_point_on_side    Index of point on side
     */
    uint register_side_point(uint patch_side_idx, uint elm_cache_map_idx, uint elem_idx, uint side_idx, uint i_point_on_side) {
        uint point_pos = i_point_on_side * n_sides() + patch_side_idx; // index of side point on patch
        int_table_(0)(point_pos) = elm_cache_map_idx;
        int_table_(1)(point_pos) = patch_side_idx;
        int_table_(2)(point_pos) = elem_idx;
        int_table_(4)(point_pos) = side_idx;

        points_map_[elm_cache_map_idx] = point_pos;
        return point_pos;
    }

    template<class ElementDomain>
    NodeAccessor<spacedim> node(unsigned int i_elm, unsigned int i_n);

    template<class ElementDomain>
    unsigned int n_mesh_entities();

    /// Set number of elements and points as permanent
    inline void make_permanent_mesh_items() {
        n_sides_.make_permanent();
        n_points_.make_permanent();
    }
//protected:

    /**
     * Hold integer values of quadrature points of defined operations.
     *
     * Table contains following rows:
     *  0: Index of quadrature point in ElementCacheMap
     *  1: Index of element (bulk PPV) / side (side PPV) in PatchOp::result_ table to which quadrature point is relevant
     *  2: Element idx in Mesh
     *   - last four rows are allocated only for side point table
     *  3: Index of side in element - short vector, size of column = number of sides
     *  4: Index of side in element - long vector, size of column = number of points
     *  5. Index of element in PatchOp::result_ table to which side is relevant
     * Number of used rows is given by n_points_.
     */
    IntTableArena int_table_;

    /// Set size and type of rows of int_table_, value is set implicitly in constructor of descendants
    std::vector<OpSizeType> int_sizes_;

    ElemsDimOata<spacedim> *elems_dim_data_;  ///< Number and list of elements on patch
    RevertibleValue n_points_;                ///< Number of points in patch
    RevertibleValue n_sides_;                 ///< Number of sides in patch
    uint i_side_;                             ///< Index of registered side in table, helper value used during patch creating
    std::vector<uint> points_map_;            ///< Map of point patch indices to PatchOp::result_ and int_table_ tables
	std::vector<Side> side_list_;             ///< List of sides on patch
};



#endif /* PATCH_POINT_VALUES_HH_ */
