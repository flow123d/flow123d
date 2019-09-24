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
#include "fields/composed_quadrature.hh"
#include "mesh/side_impl.hh"
#include "mesh/sides.h"


/******************************************************************************
 * Implementation of BulkSubQuad methods.
 */

template <unsigned int dim>
Range< BulkPointAccessor<dim> > BulkSubQuad<dim>::points() const {
	auto bgn_it = make_iter<BulkPointAccessor<dim>>( BulkPointAccessor<dim>(*this, 0) );
	auto end_it = make_iter<BulkPointAccessor<dim>>( BulkPointAccessor<dim>(*this, point_indices_.size()) );
	return Range<BulkPointAccessor<dim>>(bgn_it, end_it);
}


/******************************************************************************
 * Implementation of SideSubQuad methods.
 */

template <unsigned int dim>
Range< PointAccessor<dim> > SideSubQuad<dim>::points(const Side &side) const {
	auto permutation = side.element()->permutation_idx(side.side_idx());
	auto bgn_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(c_quad_, this->point_indices_[permutation][side.side_idx()]) );
	auto end_it = make_iter<PointAccessor<dim>>( PointAccessor<dim>(c_quad_, this->point_indices_[permutation][side.side_idx()+1]) );
	return Range<PointAccessor<dim>>(bgn_it, end_it);
}


/******************************************************************************
 * Implementation of BulkPointAccessor methods.
 */

template <unsigned int dim>
arma::vec::fixed<dim> BulkPointAccessor<dim>::loc_coords()
{
    return c_quad().local_points_[ this->point_set_idx() ];
}

template <unsigned int dim>
arma::vec3 BulkPointAccessor<dim>::coords()
{}


/******************************************************************************
 * Implementation of PointAccessor methods.
 */

template <unsigned int dim>
arma::vec::fixed<dim> PointAccessor<dim>::loc_coords()
{
    return c_quad_->local_points_[ this->idx() ];
}

template <unsigned int dim>
arma::vec3 PointAccessor<dim>::coords()
{}


/******************************************************************************
 * Explicit instantiation of templates
 */

template class BulkSubQuad<1>;
template class BulkSubQuad<2>;
template class BulkSubQuad<3>;

template class SideSubQuad<1>;
template class SideSubQuad<2>;
template class SideSubQuad<3>;

template class BulkPointAccessor<1>;
template class BulkPointAccessor<2>;
template class BulkPointAccessor<3>;

template class PointAccessor<1>;
template class PointAccessor<2>;
template class PointAccessor<3>;
