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
 * @file    field_evaluator.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_EVALUATOR_HH_
#define FIELD_EVALUATOR_HH_


#include <vector>
#include <memory>
#include <armadillo>
#include "mesh/range_wrapper.hh"
#include "system/asserts.hh"

class FieldCommon;
class Side;
template <unsigned int dim> class BulkPointSet;
template <unsigned int dim> class SidePointSet;
template <unsigned int dim> class PointAccessor;
template <unsigned int dim> class Quadrature;
template <int spacedim> class ElementAccessor;


template <unsigned int dim>
class FieldEvaluator {
public:
    /// Constructor
	FieldEvaluator();

	/// Getter for bulk_set_.
	inline BulkPointSet<dim> bulk_point_set() const {
	    ASSERT(bulk_set_.f_eval_ != nullptr).error("Uninitialized bulk point set!\n");
		return bulk_set_;
	}

	/// Getter for side_set_.
	inline SidePointSet<dim> side_point_set() const {
	    ASSERT(side_set_.f_eval_ != nullptr).error("Uninitialized side point set!\n");
		return side_set_;
	}

    /// Registers point set from quadrature and associates to it some fields.
    /// Returns an object referencing to the FieldEvaluator and list of its points.
    BulkPointSet<dim> add_bulk_fields(const Quadrature<dim> &, std::vector<FieldCommon *>);

    /// The same as add_bulk_fields but for points on side.
    SidePointSet<dim> add_side_fields(const Quadrature<dim-1> &, std::vector<FieldCommon *>);

    /// Calls reinit() for all registered fields with particular point sets.
    /// The point sets are created for each field by joining bulk and side points
    /// respecting the side permutation.
    void reinit(ElementAccessor<3> &);

    /// Returns range loop over all points
    Range<PointAccessor<dim>> points_range() const;

    /// Returns range loop over all bulk points
    Range<PointAccessor<dim>> bulk_range() const;

    /// Returns range loop over all side points
    Range<PointAccessor<dim>> sides_range() const;

    /// Returns range of points for given side and its permutation
    Range<PointAccessor<dim>> side_range(const Side &) const;

private:
    BulkPointSet<dim> bulk_set_;
    SidePointSet<dim> side_set_;

    std::vector<arma::vec::fixed<dim>> local_points_;  ///< Local coords of points vector
    std::vector<unsigned int> bulk_range_;             ///< Range of bulk points in previous vector (vector size = 2, holds 'begin' and 'end' index)
    std::vector<unsigned int> side_ranges_;            ///< Ranges of side points (vector size = dim+2, indexes for i-th side: begin = i, end = i+1)

    std::vector<FieldCommon *> bulk_field_vec_;
    std::vector<FieldCommon *> side_field_vec_;

    friend class BulkPointSet<dim>;
    friend class SidePointSet<dim>;
    friend class PointAccessor<dim>;

};


#endif /* FIELD_EVALUATOR_HH_ */
