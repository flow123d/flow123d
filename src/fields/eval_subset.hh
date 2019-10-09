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


class EvalSubset {
public:
    /// Constructor
	EvalSubset() : eval_points_(nullptr) {}

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
        return point_indices_.size();
    }

    /// Return vector of point indices of given permutation
    const std::vector<int> &get_point_indices(unsigned int permutation) const;

    /// Returns range of bulk local points for appropriate cell accessor
    Range< BulkPoint > points(const DHCellAccessor &cell) const;

    /// Temporary method for testing
    void print_side_points(unsigned int permutation);

private:
    /// Empty subset is need in default constructors of BulkPoint and SidePoint
    static const EvalSubset dummy_subset;
    /// Pointer to composed quadrature
    const EvalPoints *eval_points_;
    /// Indices into full set of local indices in the composed quadrature, for every possible permutation.
    std::vector< std::vector<int> > point_indices_;
    /// Number of sides (value 0 indicates bulk set)
    unsigned int n_sides_;

    friend class EvalPoints;
    friend class BulkPoint;
};


class BulkPoint {
public:
    /// Default constructor
	BulkPoint()
    : local_point_idx_(0) {}

    /// Constructor
	BulkPoint(DHCellAccessor dh_cell, EvalSubset bulk_subset, unsigned int loc_point_idx)
    : subset_(bulk_subset), local_point_idx_(loc_point_idx) {}

    /// Getter of composed quadrature
    inline const EvalPoints &eval_points() const {
        return subset_.eval_points();
    }

    /// Index of point within EvalPoints
    inline unsigned int point_set_idx() const {
        return subset_.point_indices_[0][local_point_idx_];
    }

    /// Local coordinates within element
    arma::vec loc_coords() const;

    // Global coordinates within element
    //arma::vec3 coords() const;

    /// Iterates to next point.
    inline void inc() {
    	local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const BulkPoint& other) {
    	return (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// DOF handler accessor of element.
    DHCellAccessor dh_cell_;
    /// Reference to bulk point set.
    EvalSubset subset_;
    /// Index of the local point in bulk point set.
    unsigned int local_point_idx_;
};


/*template <unsigned int dim>
class SidePointAccessor {
public:
    /// Default constructor
	SidePointAccessor()
    : local_point_idx_(0), permutation_(0) {}

    /// Constructor
	SidePointAccessor(EvalSubset<dim> side_points, unsigned int local_point_idx, unsigned int perm)
    : side_points_(side_points), local_point_idx_(local_point_idx), permutation_(perm) {}

    /// Getter of composed quadrature
    inline const EvalPoints<dim> &eval_points() const {
        return *side_points_.eval_points_;
    }

    // Index of point within EvalPoints
    inline unsigned int point_set_idx() const {
        return side_points_.point_indices_[permutation_][local_point_idx_];
    }

    // Local coordinates within element
    arma::vec::fixed<dim> loc_coords();

    // Global coordinates within element
    arma::vec3 coords();

    /// Iterates to next point.
    inline void inc() {
    	local_point_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const SidePointAccessor<dim>& other) {
    	return (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// Pointer to side point set
    EvalSubset<dim> side_points_;
    /// Index of the local point in the composed quadrature.
    unsigned int local_point_idx_;
    /// Permutation of nodes on sides (equivalent with \p RefElement::side_permutations)
    unsigned int permutation_;
};*/


#endif /* POINT_SETS_HH_ */
