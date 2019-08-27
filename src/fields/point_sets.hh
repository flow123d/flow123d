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
#include "mesh/range_Wrapper.hh"


class Side;
template <unsigned int dim> class FieldEvaluator;
template <unsigned int dim> class PointAccessor;


template <unsigned int dim>
class BulkPointSet {
public:
    /// Constructor
	BulkPointSet() : f_eval_(nullptr) {}

    /// Getter of field evaluator
    inline std::shared_ptr<FieldEvaluator<dim>> field_evaluator() const {
        return f_eval_;
    }

    Range<PointAccessor<dim>> points() const;

private:
    std::shared_ptr<FieldEvaluator<dim>> f_eval_;

    friend class FieldEvaluator<dim>;
};


template <unsigned int dim>
class SidePointSet {
public:
    /// Constructor
	SidePointSet() : f_eval_(nullptr) {}

    /// Getter of field evaluator
    inline std::shared_ptr<FieldEvaluator<dim>> field_evaluator() const {
        return f_eval_;
    }

    /// returns points for given side and its permutation
    Range<PointAccessor<dim>> points(const Side &) const;

private:
    std::shared_ptr<FieldEvaluator<dim>> f_eval_;

    friend class FieldEvaluator<dim>;
};


template <unsigned int dim>
class PointAccessor {
public:
    /// Default constructor
    PointAccessor()
    : f_eval_(nullptr), idx_(0) {}

    /// Constructor
    PointAccessor(std::shared_ptr<FieldEvaluator<dim>> f_eval, unsigned int idx)
    : f_eval_(f_eval), idx_(idx) {}

    /// Getter of field evaluator
    inline std::shared_ptr<FieldEvaluator<dim>> field_evaluator() const {
        return f_eval_;
    }

    // Index of point within FieldEvaluator
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
    std::shared_ptr<FieldEvaluator<dim>> f_eval_;
    unsigned int idx_;
};


#endif /* POINT_SETS_HH_ */
