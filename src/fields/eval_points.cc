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

#include "mesh/ref_element.hh"
#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "quadrature/quadrature.hh"
#include <memory>


const unsigned int EvalPoints::undefined_dim = 10;

EvalPoints::EvalPoints()
: local_points_(0, 1), dim_(EvalPoints::undefined_dim)
{
	subset_starts_.push_back(0);
}

template <unsigned int dim>
std::shared_ptr<EvalSubset> EvalPoints::add_bulk(const Quadrature &quad)
{
	check_dim(quad.dim(), dim);

	std::shared_ptr<EvalSubset> bulk_set = std::make_shared<EvalSubset>(shared_from_this() );
	this->add_local_points<dim>( quad.get_points() );
	subset_starts_.push_back( this->size() );
    return bulk_set;
}

template <unsigned int dim>
std::shared_ptr<EvalSubset> EvalPoints::add_side(const Quadrature &quad)
{
	check_dim(quad.dim()+1, dim);
	unsigned int old_data_size=this->size(), new_data_size; // interval of side subset data
	unsigned int points_per_side = quad.make_from_side<dim>(0, 0).get_points().n_vals();

	std::shared_ptr<EvalSubset> side_set = std::make_shared<EvalSubset>(shared_from_this(), RefElement<dim>::n_side_permutations, points_per_side);
	unsigned int*** perm_indices = side_set->perm_indices_;

    // permutation 0
    for (unsigned int i=0; i<dim+1; ++i) {  // sides
        Quadrature high_dim_q = quad.make_from_side<dim>(i, 0);
        this->add_local_points<dim>( high_dim_q.get_points() );
    }
    new_data_size = this->size();
    subset_starts_.push_back( new_data_size );
    unsigned int i_data=old_data_size;
    for (unsigned int i_side=0; i_side<dim+1; ++i_side) {
        for (unsigned int i_point=0; i_point<points_per_side; ++i_point) {
        	perm_indices[i_side][0][i_point] = i_data;
        	++i_data;
        }
    }

    // permutation 1...N
    for (unsigned int i_perm=1; i_perm<RefElement<dim>::n_side_permutations; ++i_perm) {
        for (unsigned int i_side=0; i_side<dim+1; ++i_side) {
            Quadrature high_dim_q = quad.make_from_side<dim>(i_side, i_perm);
            const Armor::array & quad_points = high_dim_q.get_points();
            for (unsigned int i_point=0; i_point<quad_points.n_vals(); ++i_point) {
            	perm_indices[i_side][i_perm][i_point] = this->find_permute_point<dim>( quad_points.get<dim>(i_point).arma(), old_data_size, new_data_size );
            }
        }
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

template <unsigned int dim>
unsigned int EvalPoints::find_permute_point(arma::vec coords, unsigned int data_begin, unsigned int data_end) {
	for (unsigned int loc_idx=data_begin; loc_idx<data_end; ++loc_idx) {
	    // Check if point exists in local points vector.
        if ( arma::norm(coords-local_points_.get<dim>(loc_idx).arma(), 2) < 4*std::numeric_limits<double>::epsilon() ) return loc_idx;
    }

	ASSERT(false);
    return 0;
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


template std::shared_ptr<EvalSubset> EvalPoints::add_bulk<1>(const Quadrature &);
template std::shared_ptr<EvalSubset> EvalPoints::add_bulk<2>(const Quadrature &);
template std::shared_ptr<EvalSubset> EvalPoints::add_bulk<3>(const Quadrature &);
template std::shared_ptr<EvalSubset> EvalPoints::add_side<1>(const Quadrature &);
template std::shared_ptr<EvalSubset> EvalPoints::add_side<2>(const Quadrature &);
template std::shared_ptr<EvalSubset> EvalPoints::add_side<3>(const Quadrature &);
template void EvalPoints::add_local_points<1>(const Armor::array &);
template void EvalPoints::add_local_points<2>(const Armor::array &);
template void EvalPoints::add_local_points<3>(const Armor::array &);
template unsigned int EvalPoints::find_permute_point<1>(arma::vec, unsigned int, unsigned int);
template unsigned int EvalPoints::find_permute_point<2>(arma::vec, unsigned int, unsigned int);
template unsigned int EvalPoints::find_permute_point<3>(arma::vec, unsigned int, unsigned int);
