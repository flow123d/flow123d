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
 * @file    eval_points.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef FIELD_EVALUATOR_HH_
#define FIELD_EVALUATOR_HH_


#include <vector>
#include <memory>
#include <armadillo>
#include "mesh/range_wrapper.hh"
#include "fields/eval_subset.hh"
#include "system/asserts.hh"

class Side;
template <unsigned int dim> class Quadrature;
template <int spacedim> class ElementAccessor;


/**
 * Class holds local coordinations of evaluating points (bulk and sides).
 */
class EvalPoints {
public:
    /// Constructor
	EvalPoints();

    /// Return size of evaluation points object (number of points).
    inline unsigned int size() const {
        return local_points_.size() / dim_;
    }

    /// Return local coordinates of given local point.
    inline arma::vec local_point(unsigned int local_point_idx) const {
    	ASSERT_LT_DBG(local_point_idx, this->size());
    	std::vector<double> loc_point_vec(dim_);
    	for (unsigned int i=0; i<dim_; ++i) loc_point_vec[i] = local_points_[local_point_idx*dim_ + i];
    	return arma::vec(loc_point_vec);
    }

    /// Return dimension of stored evaluate points
    inline unsigned int point_dim() const {
        return dim_;
    }

    /**
     * Registers point set from quadrature.
     * Returns an object referencing to the EvalPoints and list of its points.
     */
    template <unsigned int dim>
    EvalSubset add_bulk(const Quadrature<dim> &);

    /// The same as add_bulk but for points on sides.
    template <unsigned int dim>
	EvalSubset add_side(const Quadrature<dim-1> &);

private:
    /// Adds coords of local point if point doesn't exist in local_points_ vector, returns its index in vector.
    unsigned int add_local_point(arma::vec coords);

    std::vector<double> local_points_;  ///< Local coords of points vector
    unsigned int dim_;                  ///< Dimension of local points

};


#endif /* FIELD_EVALUATOR_HH_ */
