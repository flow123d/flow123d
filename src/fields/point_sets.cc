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
 * @file    point_sets.cc
 * @brief
 * @author  David Flanderka
 */

#include <armadillo>
#include "fields/point_sets.hh"
#include "fields/field_evaluator.hh"
#include "mesh/sides.h"


/******************************************************************************
 * Implementation of BulkPointSet methods.
 */

template <unsigned int dim>
Range< PointAccessor<dim> > BulkPointSet<dim>::points() const {
	auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(f_eval_, f_eval_->bulk_range_[0]) );
	auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(f_eval_, f_eval_->bulk_range_[1]) );
	return Range<PointAccessor<dim>>(bgn_it, end_it);
}


/******************************************************************************
 * Implementation of SidePointSet methods.
 */

template <unsigned int dim>
Range< PointAccessor<dim> > SidePointSet<dim>::points(const Side &side) const {
	auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(f_eval_, f_eval_->side_ranges_[side.side_idx()]) );
	auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(f_eval_, f_eval_->side_ranges_[side.side_idx()+1]) );
	return Range<PointAccessor<dim>>(bgn_it, end_it);
}


/******************************************************************************
 * Implementation of PointAccessor methods.
 */

template <unsigned int dim>
arma::vec::fixed<dim> PointAccessor<dim>::loc_coords()
{
    return f_eval_->local_points_[idx_];
}

template <unsigned int dim>
arma::vec3 PointAccessor<dim>::coords()
{}


/******************************************************************************
 * Explicit instantiation of templates
 */

template class BulkPointSet<1>;
template class BulkPointSet<2>;
template class BulkPointSet<3>;

template class SidePointSet<1>;
template class SidePointSet<2>;
template class SidePointSet<3>;

template class PointAccessor<1>;
template class PointAccessor<2>;
template class PointAccessor<3>;
