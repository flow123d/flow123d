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
#include "mesh/ref_element.hh"


class Side;
template <unsigned int dim> class ComposedQuadrature;
template <unsigned int dim> class BulkPointAccessor;
template <unsigned int dim> class PointAccessor;


template <unsigned int dim>
class BulkSubQuad {
public:
    /// Constructor
	BulkSubQuad() : c_quad_(nullptr) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    /// Returnx range of bulk local points.
    Range<BulkPointAccessor<dim>> points() const;

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    /// Indices into full set of local indices in the composed quadrature.
    std::vector<int> point_indices_;

    friend class ComposedQuadrature<dim>;
    friend class BulkPointAccessor<dim>;
};


template <unsigned int dim>
class SideSubQuad {
public:
    /// Constructor
	SideSubQuad() : c_quad_(nullptr), point_indices_(RefElement<dim>::n_side_permutations) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    /// Returns range of side local points for given side and its permutation.
    Range<PointAccessor<dim>> points(const Side &side) const;

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    /// Indices into full set of local indices in the composed quadrature, for every possible permuation.
    std::vector< std::vector<int> > point_indices_;

    friend class ComposedQuadrature<dim>;
};


template <unsigned int dim>
class BulkPointAccessor {
public:
    /// Default constructor
	BulkPointAccessor()
    : local_point_idx_(0) {}

    /// Constructor
	BulkPointAccessor(BulkSubQuad<dim> bulk_points, unsigned int loc_point_idx)
    : bulk_points_(bulk_points), local_point_idx_(loc_point_idx) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *bulk_points_.c_quad_;
    }

    // Index of point within ComposedQuadrature
    inline unsigned int point_set_idx() const {
        return bulk_points_.point_indices_[local_point_idx_];
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
    bool operator==(const BulkPointAccessor<dim>& other) {
    	return (local_point_idx_ == other.local_point_idx_);
    }

private:
    /// Pointer to bulk point set.
    BulkSubQuad<dim> bulk_points_;
    /// Index of the local point in bulk point set.
    unsigned int local_point_idx_;
};


template <unsigned int dim>
class PointAccessor {
public:
    /// Default constructor
    PointAccessor()
    : c_quad_(nullptr), idx_(0) {}

    /// Constructor
    PointAccessor(const ComposedQuadrature<dim> *c_quad, unsigned int idx)
    : c_quad_(c_quad), idx_(idx) {}

    /// Getter of composed quadrature
    inline const ComposedQuadrature<dim> &c_quad() const {
        return *c_quad_;
    }

    // Index of point within ComposedQuadrature
    inline unsigned int idx() const {
        return idx_;
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
    	return (idx_ == other.idx_);
    }

private:
    /// Pointer to composed quadrature
    const ComposedQuadrature<dim> *c_quad_;
    /// Index of the local point in the composed quadrature.
    unsigned int idx_;
};


#endif /* POINT_SETS_HH_ */
