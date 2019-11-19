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
: local_points_(0, 1), dim_(EvalPoints::undefined_dim)
{
	block_starts_.push_back(0);
}

template <unsigned int dim>
EvalSubset EvalPoints::add_bulk(const Quadrature &quad)
{
	check_dim(quad.dim(), dim);

	EvalSubset bulk_set(shared_from_this() );
	this->add_local_points<dim>( quad.get_points() );
	block_starts_.push_back( this->size() );
    return bulk_set;
}

template <unsigned int dim>
EvalSubset EvalPoints::add_side(const Quadrature &quad)
{
	check_dim(quad.dim()+1, dim);

	EvalSubset side_set(shared_from_this(), RefElement<dim>::n_side_permutations);

    for (unsigned int j=0; j<RefElement<dim>::n_side_permutations; ++j) { // permutations
        for (unsigned int i=0; i<dim+1; ++i) {  // sides
            Quadrature high_dim_q = quad.make_from_side<dim>(i, j);
            this->add_local_points<dim>( high_dim_q.get_points() );
        }
        block_starts_.push_back( this->size() );
    }

    return side_set;
}

template <unsigned int dim>
void EvalPoints::add_local_points(const Armor::array & quad_points) {
    unsigned int local_points_old_size = local_points_.n_vals();
    local_points_.resize(quad_points.n_vals() + local_points_old_size);
    for (unsigned int i=0; i<quad_points.n_vals(); ++i) {
        local_points_.get<dim>(i+local_points_old_size) = quad_points.get<dim>(i);
	}
}


unsigned int EvalPoints::check_dim(unsigned int quad_dim, unsigned int obj_dim) {
	ASSERT_EQ(quad_dim, obj_dim);
    if (this->dim_ == EvalPoints::undefined_dim) {
        this->dim_ = quad_dim;
        local_points_ = Armor::array(0, this->dim_);
    } else
        ASSERT_EQ(this->dim_, quad_dim);
    return this->dim_;
}


template EvalSubset EvalPoints::add_bulk<1>(const Quadrature &);
template EvalSubset EvalPoints::add_bulk<2>(const Quadrature &);
template EvalSubset EvalPoints::add_bulk<3>(const Quadrature &);
template EvalSubset EvalPoints::add_side<1>(const Quadrature &);
template EvalSubset EvalPoints::add_side<2>(const Quadrature &);
template EvalSubset EvalPoints::add_side<3>(const Quadrature &);
template void EvalPoints::add_local_points<1>(const Armor::array &);
template void EvalPoints::add_local_points<2>(const Armor::array &);
template void EvalPoints::add_local_points<3>(const Armor::array &);
