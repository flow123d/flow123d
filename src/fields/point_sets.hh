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
 * @file    point_sets.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef POINT_SETS_HH_
#define POINT_SETS_HH_

#include <memory>
#include <armadillo>
#include "mesh/range_wrapper.hh"


class Side;
template <unsigned int dim> class ComposedQuadrature;
template <unsigned int dim> class PointAccessor;


template <unsigned int dim>
class BulkSubQuad {
public:
    /// Constructor
	BulkSubQuad() : c_quad_(nullptr), point_indices_(2, 0) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    /// Returnx range of bulk local points.
    Range<PointAccessor<dim>> points() const;

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    ///< Range of bulk points in ComposedQuadrature::local_points_ vector (vector size = 2, holds 'begin' and 'end' index)
    std::vector<int> point_indices_;

    friend class ComposedQuadrature<dim>;
};


template <unsigned int dim>
class SideSubQuad {
public:
    /// Constructor
	SideSubQuad() : c_quad_(nullptr), point_indices_(dim+2, 0) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    /// Returns range of side local points for given side and its permutation.
    Range<PointAccessor<dim>> points(const Side &side, const unsigned int side_permutations[dim]) const;

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    /// Ranges of side points (vector size = dim+2, indexes for i-th side: begin = i, end = i+1)
    std::vector<int> point_indices_;

    friend class ComposedQuadrature<dim>;
};


template <unsigned int dim>
class PointAccessor {
public:
    /// Default constructor
    PointAccessor()
    : c_quad_(nullptr), idx_(0), perm_shift_(0) {}

    /// Constructor
    PointAccessor(const ComposedQuadrature<dim> *c_quad, unsigned int idx, std::array<int, dim> permutation)
    : c_quad_(c_quad), idx_(idx), permutation_(permutation), perm_shift_(idx) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    // Index of point within ComposedQuadrature
    inline unsigned int idx() const {
        return idx_ + permutation_[idx_-perm_shift_];
    }

    // Local coordinates within element
    arma::vec::fixed<dim> loc_coords();

    // Global coordinates within element
    arma::vec3 coords();

    /// Iterates to next point.
    inline void inc() {
    	idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const PointAccessor<dim>& other) {
    	return (idx_ == other.idx_) && (permutation_ == other.permutation_);
    }

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    /// Index of the local point in the composed quadrature.
    unsigned int idx_;
    /// Holds shift of value computing from permutation (used only for side point accessors)
    std::array<int, dim> permutation_;
    /// Helper data member of range method, holds shift between indexes of range and permutation array.
    unsigned int perm_shift_;
};


#endif /* POINT_SETS_HH_ */
