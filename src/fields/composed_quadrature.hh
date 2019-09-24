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
 * @file    composed_quadrature.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_EVALUATOR_HH_
#define FIELD_EVALUATOR_HH_


#include <vector>
#include <memory>
#include <armadillo>
#include "mesh/range_wrapper.hh"
#include "fields/point_sets.hh"
#include "system/asserts.hh"

class Side;
template <unsigned int dim> class Quadrature;
template <int spacedim> class ElementAccessor;


/**
 * Class holds local coordinations of evaluating points (bulk and sides).
 *
 * TODO: rename to LocalPointSet
 */
template <unsigned int dim>
class ComposedQuadrature {
public:
    /// Constructor
	ComposedQuadrature();

	/// Getter for bulk sub quadrature.
	inline BulkSubQuad<dim> bulk_quad() const {
	    ASSERT(bulk_set_.c_quad_ != nullptr).error("Uninitialized bulk point set!\n");
		return bulk_set_;
	}

	/// Getter for side sub quadrature.
	inline SideSubQuad<dim> side_quad() const {
	    ASSERT(side_set_.c_quad_ != nullptr).error("Uninitialized side point set!\n");
		return side_set_;
	}

    /// Return size of composed quadrature (number of points).
    inline unsigned int size() const {
        return local_points_.size();
    }

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the ComposedQuadrature and list of its points.
     */
	BulkSubQuad<dim> add_bulk(const Quadrature<dim> &);

    /// The same as add_bulk but for points on side.
	SideSubQuad<dim> add_side(const Quadrature<dim-1> &);

    /// Returns range loop over all bulk points
    Range<BulkPointAccessor<dim>> bulk_range() const;

    /// Returns range loop over all side points
    //Range<PointAccessor<dim>> sides_range() const;

    /// Returns range of points for given side and its permutation
    Range<SidePointAccessor<dim>> side_range(const Side &side) const;

private:
    /// Adds coords of local point if point doesn't exist in local_points_ vector, returns its index in vector.
    unsigned int add_local_point(arma::vec::fixed<dim> coords);

    BulkSubQuad<dim> bulk_set_;  ///< Handler to bulk local points.
    SideSubQuad<dim> side_set_;  ///< Handler to sides local points.

    std::vector<arma::vec::fixed<dim>> local_points_;  ///< Local coords of points vector

    friend class BulkSubQuad<dim>;
    friend class SideSubQuad<dim>;
    friend class BulkPointAccessor<dim>;
    friend class SidePointAccessor<dim>;

};


#endif /* FIELD_EVALUATOR_HH_ */
