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


const unsigned int EvalPoints::undefined_dim = 10;

EvalPoints::EvalPoints()
: dim_eval_points_({DimEvalPoints(1), DimEvalPoints(2), DimEvalPoints(3)})
{}

template <unsigned int dim>
std::shared_ptr<BulkIntegral> EvalPoints::add_bulk(const Quadrature &quad)
{
    ASSERT_EQ(dim, quad.dim());
    std::shared_ptr<BulkIntegral> bulk_integral = std::make_shared<BulkIntegral>(shared_from_this(), dim);
    dim_eval_points_[dim-1].add_local_points<dim>( quad.get_points() );
    dim_eval_points_[dim-1].add_subset();
    return bulk_integral;
}

template <unsigned int dim>
std::shared_ptr<EdgeIntegral> EvalPoints::add_edge(const Quadrature &quad)
{
    ASSERT_EQ(dim, quad.dim()+1);
	unsigned int old_data_size=this->size(dim), new_data_size; // interval of side subset data
	unsigned int points_per_side = quad.make_from_side<dim>(0, 0).get_points().size();
	unsigned int n_side_permutations = RefElement<dim>::n_side_permutations;

	std::shared_ptr<EdgeIntegral> edge_integral = std::make_shared<EdgeIntegral>(shared_from_this(), dim, n_side_permutations, points_per_side);
	unsigned int*** perm_indices = edge_integral->perm_indices_;

    // permutation 0
    for (unsigned int i=0; i<dim+1; ++i) {  // sides
        Quadrature high_dim_q = quad.make_from_side<dim>(i, 0);
        dim_eval_points_[dim-1].add_local_points<dim>( high_dim_q.get_points() );
    }
    dim_eval_points_[dim-1].add_subset();
    new_data_size = this->size(dim);
    unsigned int i_data=old_data_size;
    for (unsigned int i_side=0; i_side<dim+1; ++i_side) {
        for (unsigned int i_point=0; i_point<points_per_side; ++i_point) {
        	perm_indices[i_side][0][i_point] = i_data;
        	++i_data;
        }
    }

    // permutation 1...N
    for (unsigned int i_perm=1; i_perm<n_side_permutations; ++i_perm) {
        for (unsigned int i_side=0; i_side<dim+1; ++i_side) {
            Quadrature high_dim_q = quad.make_from_side<dim>(i_side, i_perm);
            const Armor::Array<double> & quad_points = high_dim_q.get_points();
            for (unsigned int i_point=0; i_point<quad_points.size(); ++i_point) {
                perm_indices[i_side][i_perm][i_point] = dim_eval_points_[dim-1].find_permute_point<dim>( quad_points.vec<dim>(i_point), old_data_size, new_data_size );
            }
        }
    }

    return edge_integral;
}

template <unsigned int dim>
std::shared_ptr<CouplingIntegral> EvalPoints::add_coupling(const Quadrature &quad) {
    ASSERT_EQ(dim, quad.dim()+1);

    std::shared_ptr<BulkIntegral> bulk_integral = this->add_bulk<dim-1>(quad);
    std::shared_ptr<EdgeIntegral> edge_integral = this->add_edge<dim>(quad);
    return std::make_shared<CouplingIntegral>(edge_integral, bulk_integral);
}

template <unsigned int dim>
std::shared_ptr<BoundaryIntegral> EvalPoints::add_boundary(const Quadrature &quad) {
    ASSERT_EQ(dim, quad.dim()+1);
    return std::make_shared<BoundaryIntegral>(this->add_edge<dim>(quad));
}

EvalPoints::DimEvalPoints::DimEvalPoints(unsigned int dim)
: local_points_(dim), n_subsets_(0), dim_(dim)
{
	subset_starts_[0] = 0;
	local_points_.reinit(EvalPoints::max_subsets*EvalPoints::max_subset_points);
}

template <unsigned int dim>
void EvalPoints::DimEvalPoints::add_local_points(const Armor::Array<double> & quad_points) {
    unsigned int local_points_old_size = local_points_.size();
    local_points_.resize(quad_points.size() + local_points_old_size);
    for (unsigned int i=0; i<quad_points.size(); ++i) {
        local_points_.set(i+local_points_old_size) = quad_points.vec<dim>(i);
	}
}

template <unsigned int dim>
unsigned int EvalPoints::DimEvalPoints::find_permute_point(arma::vec coords, unsigned int data_begin, unsigned int data_end) {
	for (unsigned int loc_idx=data_begin; loc_idx<data_end; ++loc_idx) {
	    // Check if point exists in local points vector.
        if ( arma::norm(coords-local_points_.vec<dim>(loc_idx), 2) < 4*std::numeric_limits<double>::epsilon() ) return loc_idx;
    }

	ASSERT(false);
    return 0;
}

void EvalPoints::DimEvalPoints::add_subset() {
    ASSERT_LT_DBG(n_subsets_, EvalPoints::max_subsets).error("Maximal number of subsets exceeded!\n");

    n_subsets_++;
    subset_starts_[n_subsets_] = this->size();
}


template std::shared_ptr<BulkIntegral> EvalPoints::add_bulk<1>(const Quadrature &);
template std::shared_ptr<BulkIntegral> EvalPoints::add_bulk<2>(const Quadrature &);
template std::shared_ptr<BulkIntegral> EvalPoints::add_bulk<3>(const Quadrature &);
template std::shared_ptr<EdgeIntegral> EvalPoints::add_edge<1>(const Quadrature &);
template std::shared_ptr<EdgeIntegral> EvalPoints::add_edge<2>(const Quadrature &);
template std::shared_ptr<EdgeIntegral> EvalPoints::add_edge<3>(const Quadrature &);
template std::shared_ptr<CouplingIntegral> EvalPoints::add_coupling<2>(const Quadrature &);
template std::shared_ptr<CouplingIntegral> EvalPoints::add_coupling<3>(const Quadrature &);
template std::shared_ptr<BoundaryIntegral> EvalPoints::add_boundary<1>(const Quadrature &);
template std::shared_ptr<BoundaryIntegral> EvalPoints::add_boundary<2>(const Quadrature &);
template std::shared_ptr<BoundaryIntegral> EvalPoints::add_boundary<3>(const Quadrature &);
template void EvalPoints::DimEvalPoints::add_local_points<1>(const Armor::Array<double> &);
template void EvalPoints::DimEvalPoints::add_local_points<2>(const Armor::Array<double> &);
template void EvalPoints::DimEvalPoints::add_local_points<3>(const Armor::Array<double> &);
template unsigned int EvalPoints::DimEvalPoints::find_permute_point<1>(arma::vec, unsigned int, unsigned int);
template unsigned int EvalPoints::DimEvalPoints::find_permute_point<2>(arma::vec, unsigned int, unsigned int);
template unsigned int EvalPoints::DimEvalPoints::find_permute_point<3>(arma::vec, unsigned int, unsigned int);
