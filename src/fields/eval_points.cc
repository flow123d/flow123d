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
 * @file    eval_points.cc
 * @brief
 * @author  David Flanderka
 */

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "quadrature/quadrature.hh"
#include "mesh/ref_element.hh"
#include <memory>


EvalPoints::EvalPoints()
: dim_(EvalPoints::undefined_dim)
{
    block_indices_.push_back(0);
}

template <unsigned int dim>
EvalSubset EvalPoints::add_bulk(const Quadrature &quad)
{
	check_dim(quad.dim(), dim);

	EvalSubset bulk_set(this);

    const Armor::array & quad_points = quad.get_points();
    for (uint i=0; i<quad_points.n_vals(); ++i)
        this->add_local_point(quad_points.get<dim>(i).arma());
    block_indices_.push_back( this->size() );
    return bulk_set;
}

template <unsigned int dim>
EvalSubset EvalPoints::add_side(const Quadrature &quad)
{
	check_dim(quad.dim()+1, dim);

	EvalSubset side_set(this, RefElement<dim>::n_side_permutations);

    for (unsigned int j=0; j<RefElement<dim>::n_side_permutations; ++j) { // permutations
        for (unsigned int i=0; i<dim+1; ++i) {  // sides
            Quadrature high_dim_q = quad.make_from_side<dim>(i, j);
            const Armor::array & quad_points = high_dim_q.get_points();
            for (uint k=0; k<quad_points.n_vals(); ++k) {
            	this->add_local_point(quad_points.get<dim>(k).arma());
            }
        }
        block_indices_.push_back( this->size() );
    }

    return side_set;
}

unsigned int EvalPoints::add_local_point(arma::vec coords) {
	// Add new point
	for (unsigned int i=0; i<dim_; ++i) local_points_.push_back(coords[i]);
    return this->size()-1;
}


unsigned int EvalPoints::check_dim(unsigned int quad_dim, unsigned int obj_dim) {
	ASSERT_EQ(quad_dim, obj_dim);
    if (this->dim_ == EvalPoints::undefined_dim)
        this->dim_ = quad_dim;
    else
        ASSERT_EQ(this->dim_, quad_dim);
    return this->dim_;
}


template EvalSubset EvalPoints::add_bulk<1>(const Quadrature &);
template EvalSubset EvalPoints::add_bulk<2>(const Quadrature &);
template EvalSubset EvalPoints::add_bulk<3>(const Quadrature &);
template EvalSubset EvalPoints::add_side<1>(const Quadrature &);
template EvalSubset EvalPoints::add_side<2>(const Quadrature &);
template EvalSubset EvalPoints::add_side<3>(const Quadrature &);
