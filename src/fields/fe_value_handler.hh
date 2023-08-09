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
 * @file    fe_value_handler.hh
 * @brief
 */

#ifndef FE_VALUE_HANDLER_HH_
#define FE_VALUE_HANDLER_HH_

#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "mesh/accessors.hh"


/**
 * Compute real coordinates and weights (use QGauss) for given element.
 */
template <int elemdim, int spacedim>
unsigned int compute_fe_quadrature(std::vector<arma::vec::fixed<3>> & q_points, std::vector<double> & q_weights,
		const ElementAccessor<spacedim> &elm, unsigned int order=3)
{
	static const double weight_coefs[] = { 1., 1., 2., 6. };

	QGauss qgauss(elemdim, order);
	arma::mat map_mat = MappingP1<elemdim,spacedim>::element_map(elm);

	for(unsigned i=0; i<qgauss.size(); ++i) {
		q_weights[i] = qgauss.weight(i)*weight_coefs[elemdim];
		q_points[i] = MappingP1<elemdim,spacedim>::project_unit_to_real(RefElement<elemdim>::local_to_bary(qgauss.point<elemdim>(i)), map_mat);
	}

	return qgauss.size();
}



#endif /* FE_VALUE_HANDLER_HH_ */
