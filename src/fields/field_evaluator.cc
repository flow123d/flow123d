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
 * @file    field_evaluator.cc
 * @brief
 * @author  David Flanderka
 */

#include "fields/field_evaluator.hh"
#include "fields/point_sets.hh"
#include "quadrature/quadrature.hh"
#include <memory>


template <unsigned int dim>
FieldEvaluator<dim>::FieldEvaluator()
: bulk_range_(2, 0), side_ranges_(dim+2, 0) {}

template <unsigned int dim>
BulkPointSet<dim> FieldEvaluator<dim>::add_bulk_fields(const Quadrature<dim> &quad, std::vector<FieldCommon *> field_vec)
{
    ASSERT(bulk_set_.f_eval_==nullptr).error("Multiple initialization of bulk point set!\n");

    bulk_set_.f_eval_ = std::shared_ptr< FieldEvaluator<dim> >(this); //std::enable_shared_from_this< FieldEvaluator<dim> >::shared_from_this();
    bulk_field_vec_ = field_vec;

    bulk_range_[0] = local_points_.size();
    for (auto p : quad.get_points()) local_points_.push_back(p);
    bulk_range_[1] = local_points_.size();

    return bulk_set_;
}

template <unsigned int dim>
SidePointSet<dim> FieldEvaluator<dim>::add_side_fields(const Quadrature<dim-1> &quad, std::vector<FieldCommon *> field_vec)
{
    ASSERT(side_set_.f_eval_==nullptr).error("Multiple initialization of side point set!\n");

    side_set_.f_eval_ = std::shared_ptr< FieldEvaluator<dim> >(this); //std::enable_shared_from_this< FieldEvaluator<dim> >::shared_from_this();
    side_field_vec_ = field_vec;

    for (unsigned int i=0; i<dim+1; ++i) {
        Quadrature<dim> high_dim_q(quad, i, 0); // set correct permutation id
        side_ranges_[i] = local_points_.size();
        for (auto p : high_dim_q.get_points()) local_points_.push_back(p);
    }
    side_ranges_[dim+1] = local_points_.size();

    return side_set_;
}

template <unsigned int dim>
void FieldEvaluator<dim>::reinit(ElementAccessor<3> &elm)
{
    // Not implemented yet!
}

template <unsigned int dim>
Range< PointAccessor<dim> > FieldEvaluator<dim>::points_range() const {
    auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->bulk_point_set().field_evaluator(), 0) );
    auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->bulk_point_set().field_evaluator(), local_points_.size()) );
    return Range<PointAccessor<dim>>(bgn_it, end_it);
}

template <unsigned int dim>
Range< PointAccessor<dim> > FieldEvaluator<dim>::bulk_range() const {
    return this->bulk_point_set().points();
}

template <unsigned int dim>
Range< PointAccessor<dim> > FieldEvaluator<dim>::sides_range() const {
    auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->side_point_set().field_evaluator(), side_ranges_[0]) );
    auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(this->side_point_set().field_evaluator(), side_ranges_[dim+2]) );
    return Range<PointAccessor<dim>>(bgn_it, end_it);
}

template <unsigned int dim>
Range< PointAccessor<dim> > FieldEvaluator<dim>::side_range(const Side &side) const {
    return this->side_point_set().points(side);
}


template class FieldEvaluator<1>;
template class FieldEvaluator<2>;
template class FieldEvaluator<3>;
