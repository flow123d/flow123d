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
#include <memory>


template <unsigned int dim>
ComposedQuadrature<dim>::ComposedQuadrature()
{}

template <unsigned int dim>
BulkSubQuad<dim> ComposedQuadrature<dim>::add_bulk_quad(const Quadrature<dim> &quad)
{
    ASSERT(bulk_set_.c_quad_==nullptr).error("Multiple initialization of bulk point set!\n");

    bulk_set_.c_quad_ = this;

    bulk_set_.point_indices_[0] = local_points_.size();
    for (auto p : quad.get_points()) local_points_.push_back(p);
    bulk_set_.point_indices_[1] = local_points_.size();

    return bulk_set_;
}

template <unsigned int dim>
SideSubQuad<dim> ComposedQuadrature<dim>::add_side_quad(const Quadrature<dim-1> &quad)
{
    ASSERT(side_set_.c_quad_==nullptr).error("Multiple initialization of side point set!\n");

    side_set_.c_quad_ = this;

    for (unsigned int i=0; i<dim+1; ++i) {
        Quadrature<dim> high_dim_q(quad, i, 0); // set correct permutation id
        side_set_.point_indices_[i] = local_points_.size();
        for (auto p : high_dim_q.get_points()) local_points_.push_back(p);
    }
    side_set_.point_indices_[dim+1] = local_points_.size();

    return side_set_;
}

template <unsigned int dim>
Range< PointAccessor<dim> > ComposedQuadrature<dim>::bulk_range() const {
    return this->bulk_quad().points();
}

/*template <unsigned int dim>
Range< PointAccessor<dim> > ComposedQuadrature<dim>::sides_range() const {
    auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->side_quad().c_quad_, side_set_.point_indices_[0]) );
    auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->side_quad().c_quad_, side_set_.point_indices_[dim+1]) );
    return Range<PointAccessor<dim>>(bgn_it, end_it);
}*/

template <unsigned int dim>
Range< PointAccessor<dim> > ComposedQuadrature<dim>::side_range(const Side &side, const unsigned int side_permutations[dim]) const {
    return this->side_quad().points(side, side_permutations);
}


template class ComposedQuadrature<1>;
template class ComposedQuadrature<2>;
template class ComposedQuadrature<3>;
