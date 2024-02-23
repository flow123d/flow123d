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
 * @file    patch_point_values.cc
 * @brief   Store finite element data on the actual patch
 *          such as shape function values, gradients, Jacobian
 *          of the mapping from the reference cell etc.
 * @author  David Flanderka
 */

#include "fem/patch_point_values.hh"


template<unsigned int spacedim>
PatchPointValues<spacedim>::PatchPointValues(uint dim, PointType point_type)
: dim_(dim), point_type_(point_type), n_columns_(0) {}


template<unsigned int spacedim>
void PatchPointValues<spacedim>::initialize(uint int_cols) {
    this->reset();

	point_vals_.resize(n_columns_);
	int_vals_.resize(int_cols);
	el_vals_.resize( (dim_+1) * spacedim );
}


template<unsigned int spacedim>
void OpJacBulk<spacedim>::reinit_data() {
    // compute data
}

template<unsigned int spacedim>
void OpJacDetBulk<spacedim>::reinit_data() {
    // compute data
}

template class PatchPointValues<3>;
