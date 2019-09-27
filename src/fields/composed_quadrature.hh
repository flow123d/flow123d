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
 */
template <unsigned int dim>
class EvalPoints {
public:
    /// Constructor
	EvalPoints();

    /// Return size of composed quadrature (number of points).
    inline unsigned int size() const {
        return local_points_.size();
    }

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the EvalPoints and list of its points.
     */
    EvalSubset<dim> add_bulk(const Quadrature<dim> &);

    /// The same as add_bulk but for points on sides.
	EvalSubset<dim> add_side(const Quadrature<dim-1> &);

private:
    /// Adds coords of local point if point doesn't exist in local_points_ vector, returns its index in vector.
    unsigned int add_local_point(arma::vec::fixed<dim> coords);

    EvalSubset<dim> bulk_set_;  ///< Handler to bulk local points.
    EvalSubset<dim> side_set_;  ///< Handler to sides local points.

    std::vector<arma::vec::fixed<dim>> local_points_;  ///< Local coords of points vector

    friend class EvalSubset<dim>;

};


#endif /* FIELD_EVALUATOR_HH_ */
