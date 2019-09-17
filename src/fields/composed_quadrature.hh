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
#include "system/asserts.hh"

class Side;
template <unsigned int dim> class BulkSubQuad;
template <unsigned int dim> class SideSubQuad;
template <unsigned int dim> class PointAccessor;
template <unsigned int dim> class Quadrature;
template <int spacedim> class ElementAccessor;


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

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the ComposedQuadrature and list of its points.
     */
	BulkSubQuad<dim> add_bulk_quad(const Quadrature<dim> &);

    /// The same as add_bulk_quad but for points on side.
	SideSubQuad<dim> add_side_quad(const Quadrature<dim-1> &);

    /// Returns range loop over all bulk points
    Range<PointAccessor<dim>> bulk_range() const;

    /// Returns range loop over all side points
    Range<PointAccessor<dim>> sides_range() const;

    /// Returns range of points for given side and its permutation
    Range<PointAccessor<dim>> side_range(const Side &) const;

private:
    BulkSubQuad<dim> bulk_set_;
    SideSubQuad<dim> side_set_;

    std::vector<arma::vec::fixed<dim>> local_points_;  ///< Local coords of points vector

    friend class BulkSubQuad<dim>;
    friend class SideSubQuad<dim>;
    friend class PointAccessor<dim>;

};


#endif /* FIELD_EVALUATOR_HH_ */
