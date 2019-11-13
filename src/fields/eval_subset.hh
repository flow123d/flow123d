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
 * @file    eval_subset.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef POINT_SETS_HH_
#define POINT_SETS_HH_

#include <memory>
#include <armadillo>
#include "mesh/range_wrapper.hh"
#include "fem/dh_cell_accessor.hh"


class Side;
class EvalPoints;
class BulkPoint;
class SidePoint;


class EvalSubset {
public:
	TYPEDEF_ERR_INFO(EI_ElementIdx, unsigned int);
	DECLARE_EXCEPTION(ExcElementNotInCache,
	        << "Element of Idx: " << EI_ElementIdx::val << " is not stored in 'Field value data cache'.\n"
			   << "Value can't be computed.\n");

    /// Constructor
	EvalSubset() : eval_points_(nullptr) {}

    /// Constructor of bulk subset
	EvalSubset(const EvalPoints *eval_points);

    /// Constructor of side subset
	EvalSubset(const EvalPoints *eval_points, unsigned int n_permutations);

    /// Getter of composed quadrature
    inline const EvalPoints &eval_points() const {
        return *eval_points_;
    }

    /// Getter of n_sides
    inline const unsigned int n_sides() const {
        return n_sides_;
    }

    /// Return number of permutations
    inline const unsigned int n_permutations() const {
        return block_indices_.size();
    }

    /// Return vector of point indices of given permutation
    inline int get_block_idx(unsigned int block_idx) const {
        ASSERT_LT( block_idx, block_indices_.size() );
        return block_indices_[block_idx];
    }

    /// Returns range of bulk local points for appropriate cell accessor
    Range< BulkPoint > points(const DHCellAccessor &cell) const;

    /// Returns range of side local points for appropriate cell side accessor
    Range< SidePoint > points(const DHCellSide &cell_side) const;

private:
    /// Empty subset is need in default constructors of BulkPoint and SidePoint
    static const EvalSubset dummy_subset;
    /// Pointer to composed quadrature
    const EvalPoints *eval_points_;
    /// Indices of data blocks in EvalPoints object for all permutations.
    std::vector<int> block_indices_;
    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;
};


class BulkPoint {
public:
    /// Default constructor
	BulkPoint()
    : local_point_idx_(0) {}

    /// Constructor
	BulkPoint(DHCellAccessor dh_cell, EvalSubset bulk_subset, unsigned int loc_point_idx)
    : dh_cell_(dh_cell), subset_(bulk_subset), local_point_idx_(loc_point_idx) {}

    /// Getter of composed quadrature
    inline const EvalPoints &eval_points() const {
        return subset_.eval_points();
    }

    /// Local coordinates within element
    arma::vec loc_coords() const;

    // Global coordinates within element
    //arma::vec3 coords() const;

    /// Return index of element in data cache.
    inline unsigned int element_cache_index() const {
        return dh_cell_.element_cache_index();
    }

    /// Iterates to next point.
    inline void inc() {
    	local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
    	return (dh_cell_ == other.dh_cell_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// DOF handler accessor of element.
    DHCellAccessor dh_cell_;
    /// Reference to bulk point set.
    EvalSubset subset_;
    /// Index of the local point in bulk point set.
    unsigned int local_point_idx_;
};


class SidePoint {
public:
    /// Default constructor
	SidePoint()
    : local_point_idx_(0) {}

    /// Constructor
	SidePoint(DHCellSide cell_side, EvalSubset subset, unsigned int local_point_idx)
    : cell_side_(cell_side), subset_(subset), local_point_idx_(local_point_idx),
	  permutation_idx_( cell_side.element()->permutation_idx( cell_side_.side_idx() ) ) {}

    /// Getter of evaluation points
    inline const EvalPoints &eval_points() const {
        return subset_.eval_points();
    }

    // Local coordinates within element
    arma::vec loc_coords() const;

    // Global coordinates within element
    //arma::vec3 coords() const;

    /// Return index of element in data cache.
    inline unsigned int element_cache_index() const {
        return cell_side_.cell().element_cache_index();
    }

    // Index of permutation
    inline unsigned int permutation_idx() const {
        return permutation_idx_;
    }

    /// Iterates to next point.
    inline void inc() {
    	local_point_idx_++;
    }

    /// Return corresponds SidePoints of neighbour side of same dimension (computing of side integrals).
    SidePoint permute(DHCellSide edg_side) const;

    /// Comparison of accessors.
    bool operator==(const SidePoint& other) {
    	return (cell_side_ == other.cell_side_) && (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// DOF handler accessor of element side.
    DHCellSide cell_side_;
    /// Reference to side point set
    EvalSubset subset_;
    /// Index of the local point in the composed quadrature.
    unsigned int local_point_idx_;
    /// Permutation index corresponding with DHCellSide
    unsigned int permutation_idx_;
};


#endif /* POINT_SETS_HH_ */
