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
PatchPointValues<spacedim>::PatchPointValues(uint dim)
: dim_(dim), n_columns_(0), elements_map_(300, 0), points_map_(300, 0) {}


template<unsigned int spacedim>
void PatchPointValues<spacedim>::initialize(uint int_cols) {
    this->reset();

	point_vals_.resize(n_columns_);
	int_vals_.resize(int_cols);
	el_vals_.resize(n_columns_);
}

template<unsigned int spacedim>
ElOp<spacedim> *PatchPointValues<spacedim>::add_accessor(ElOp<spacedim> *op_accessor) {
	operation_columns_.push_back(op_accessor);
	return op_accessor;
}

namespace FeBulk {

/** Implementation of PatchPointValues methods **/
template<unsigned int spacedim>
PatchPointValues<spacedim>::PatchPointValues(uint dim)
: ::PatchPointValues<spacedim>(dim) {
    // add instances of ElOp descendants to operation_columns_ vector
    this->add_accessor( new OpCoords(this->dim_) );
    ElOp<spacedim> *el_coords_bulk = this->add_accessor( new OpElCoords(this->dim_) );
    ElOp<spacedim> *jac_bulk = this->add_accessor( new OpJac(this->dim_, el_coords_bulk) );
    this->add_accessor( new OpJacDet(this->dim_, jac_bulk) );
}


/** Implementation of OpCoords methods **/
template<unsigned int spacedim>
void OpCoords<spacedim>::reinit_data() {
    // compute data
}


/** Implementation of OpElCoords methods **/
template<unsigned int spacedim>
void OpElCoords<spacedim>::reinit_data() {
    // compute data
}


/** Implementation of OpJac methods **/
template<unsigned int spacedim>
void OpJac<spacedim>::reinit_data() {
    // compute data
}

template<unsigned int spacedim>
void OpJac<spacedim>::check_op_dependency(::PatchPointValues<spacedim> &point_vals) {
    coords_operator_->register_columns(point_vals);
}


/** Implementation of OpJacDet methods **/
template<unsigned int spacedim>
void OpJacDet<spacedim>::reinit_data() {
    // compute data
}

template<unsigned int spacedim>
void OpJacDet<spacedim>::check_op_dependency(::PatchPointValues<spacedim> &point_vals) {
    jac_operator_->register_columns(point_vals);
}

} // closing namespace FeBulk


namespace FeSide {

template<unsigned int spacedim>
PatchPointValues<spacedim>::PatchPointValues(uint dim)
: ::PatchPointValues<spacedim>(dim) {
    // add instances of ElOp descendants to operation_columns_ vector
}

} // closing namespace FeSide

template class PatchPointValues<3>;
template class ElOp<3>;
template class FeBulk::PatchPointValues<3>;
template class FeBulk::OpJac<3>;
template class FeBulk::OpJacDet<3>;
template class FeSide::PatchPointValues<3>;
