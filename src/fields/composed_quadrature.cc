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
 * @file    composed_quadrature.cc
 * @brief
 * @author  David Flanderka
 */

#include "fields/composed_quadrature.hh"
#include "fields/point_sets.hh"
#include "quadrature/quadrature.hh"
#include "mesh/ref_element.hh"
#include <memory>


EvalPoints::EvalPoints()
{}

template <unsigned int dim>
EvalSubset EvalPoints::add_bulk(const Quadrature<dim> &quad)
{
	EvalSubset bulk_set;

	this->dim_ = dim;
    bulk_set.eval_points_ = this;
    bulk_set.dim = dim;
    bulk_set.point_indices_.resize(1);

    for (auto p : quad.get_points()) bulk_set.point_indices_[0].push_back( this->add_local_point(p) );
    return bulk_set;
}

template <unsigned int dim>
EvalSubset EvalPoints::add_side(const Quadrature<dim-1> &quad)
{
	EvalSubset side_set;

	this->dim_ = dim;
    side_set.eval_points_ = this;
    side_set.dim = dim;
    side_set.point_indices_.resize(RefElement<dim>::n_side_permutations);

    for (unsigned int j=0; j<RefElement<dim>::n_side_permutations; ++j) { // permutations
        for (unsigned int i=0; i<dim+1; ++i) {  // sides
            Quadrature<dim> high_dim_q(quad, i, j);
            for (auto p : high_dim_q.get_points()) {
            	side_set.point_indices_[j].push_back( this->add_local_point(p) );
            }
        }
    }

    return side_set;
}

unsigned int EvalPoints::add_local_point(arma::vec coords) {
    // Check if point exists in local points vector.
	std::vector<double> loc_point_vec(dim_);
	for (unsigned int loc_idx=0; loc_idx<this->size(); ++loc_idx) {
		for (unsigned int i=0; i<dim_; ++i) loc_point_vec[i] = local_points_[loc_idx*dim_ + i];
        if ( arma::norm(coords-arma::vec(loc_point_vec), 2) < 4*std::numeric_limits<double>::epsilon() ) return loc_idx;
    }
	// Add new point if doesn't exist
	for (unsigned int i=0; i<dim_; ++i) local_points_.push_back(coords[i]);
    return this->size()-1;
}


template EvalSubset EvalPoints::add_bulk<1>(const Quadrature<1> &);;
template EvalSubset EvalPoints::add_bulk<2>(const Quadrature<2> &);;
template EvalSubset EvalPoints::add_bulk<3>(const Quadrature<3> &);;
template EvalSubset EvalPoints::add_side<1>(const Quadrature<0> &);;
template EvalSubset EvalPoints::add_side<2>(const Quadrature<1> &);;
template EvalSubset EvalPoints::add_side<3>(const Quadrature<2> &);;
