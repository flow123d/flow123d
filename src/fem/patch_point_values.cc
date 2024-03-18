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
void PatchPointValues<spacedim>::initialize(Quadrature &quad, uint int_cols) {
    this->reset();

	point_vals_.resize(n_columns_);
	int_vals_.resize(int_cols);
	el_vals_.resize(n_columns_);

	ref_elm_values_ = std::make_shared<RefElementValues<spacedim> >(quad, dim_);
	ref_elm_values_ ->ref_initialize(quad, dim_);
}

template<unsigned int spacedim>
ElOp<spacedim> &PatchPointValues<spacedim>::add_accessor(ElOp<spacedim> op_accessor) {
	operations_.push_back(op_accessor);
	return operations_[operations_.size()-1];
}

template<unsigned int spacedim>
uint PatchPointValues<spacedim>::register_element(arma::mat coords, uint element_patch_idx) {
    uint res_column = operations_[FeBulk::BulkOps::opElCoords].result_col();
    for (uint i_col=0; i_col<coords.n_cols; ++i_col)
        for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
            el_vals_(res_column)(n_elems_) = coords(i_row, i_col);
            ++res_column;
        }

    elements_map_[element_patch_idx] = n_elems_;
    return n_elems_++;
}

template<unsigned int spacedim>
uint PatchPointValues<spacedim>::register_side(arma::mat coords) {
    uint res_column = operations_[FeSide::SideOps::opElCoords].result_col();
    for (uint i_col=0; i_col<coords.n_cols; ++i_col)
        for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
            el_vals_(res_column)(n_elems_) = coords(i_row, i_col);
            ++res_column;
        }

    return n_elems_++;
}

template<unsigned int spacedim>
uint PatchPointValues<spacedim>::register_bulk_point(uint elem_table_row, uint value_patch_idx, uint elem_idx) {
    int_vals_(0)(n_points_) = value_patch_idx;
    int_vals_(1)(n_points_) = elem_table_row;
    int_vals_(2)(n_points_) = elem_idx;

    points_map_[value_patch_idx] = n_points_;
    return n_points_++;
}

template<unsigned int spacedim>
uint PatchPointValues<spacedim>::register_side_point(uint elem_table_row, uint value_patch_idx, uint elem_idx, uint side_idx) {
    int_vals_(0)(n_points_) = value_patch_idx;
    int_vals_(1)(n_points_) = elem_table_row;
    int_vals_(2)(n_points_) = elem_idx;
    int_vals_(3)(n_points_) = side_idx;

    points_map_[value_patch_idx] = n_points_;
    return n_points_++;
}


namespace FeBulk {

/** Implementation of PatchPointValues methods **/
template<unsigned int spacedim>
PatchPointValues<spacedim>::PatchPointValues(uint dim)
: ::PatchPointValues<spacedim>(dim) {
    // add instances of ElOp descendants to operations_ vector
    ElOp<spacedim> &coords_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {spacedim}, this->n_columns_, false,
    				nullptr, &bulk_reinit::ptop_coords) );

    this->n_columns_ += coords_bulk.n_comp();
    ElOp<spacedim> &el_coords_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {spacedim, this->dim_+1}, this->n_columns_, true) );

    this->n_columns_ += el_coords_bulk.n_comp();
    ElOp<spacedim> &jac_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {spacedim, this->dim_}, this->n_columns_, true, &bulk_reinit::elop_jac,
    				nullptr, &el_coords_bulk) );

    this->n_columns_ += jac_bulk.n_comp();
    ElOp<spacedim> &jac_det_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {1}, this->n_columns_, true, &bulk_reinit::elop_jac_det,
    				nullptr, &jac_bulk) );

    this->n_columns_ += jac_det_bulk.n_comp();
    ElOp<spacedim> &weights_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {1}, this->n_columns_, false,
    				nullptr, &bulk_reinit::ptop_weights) );

    this->n_columns_ += weights_bulk.n_comp();
    ElOp<spacedim> &JxW_bulk = this->add_accessor(
    		ElOp<spacedim>(this->dim_, {1}, this->n_columns_, false,
    				nullptr, &bulk_reinit::ptop_JxW, &weights_bulk) );
    this->n_columns_ += JxW_bulk.n_comp();
}

} // closing namespace FeBulk


namespace FeSide {

template<unsigned int spacedim>
PatchPointValues<spacedim>::PatchPointValues(uint dim)
: ::PatchPointValues<spacedim>(dim) {
    // add instances of ElOp descendants to operations_ vector
    ElOp<spacedim> &coords_side = this->add_accessor( ElOp<spacedim>(this->dim_, {spacedim}, this->n_columns_, false, nullptr, &side_reinit::ptop_coords) );
    this->n_columns_ += coords_side.n_comp();
    ElOp<spacedim> &el_coords_side = this->add_accessor( ElOp<spacedim>(this->dim_, {spacedim, this->dim_+1}, this->n_columns_, true) );
    this->n_columns_ += el_coords_side.n_comp();
    ElOp<spacedim> &jac_side = this->add_accessor( ElOp<spacedim>(this->dim_, {spacedim, this->dim_}, this->n_columns_, true, &side_reinit::elop_jac,
            nullptr,  &el_coords_side) );
    this->n_columns_ += jac_side.n_comp();
    ElOp<spacedim> &jac_det_side = this->add_accessor( ElOp<spacedim>(this->dim_, {1}, this->n_columns_, true, &side_reinit::elop_jac_det, nullptr, &jac_side) );
    this->n_columns_ += jac_det_side.n_comp();
    ElOp<spacedim> &weights_side = this->add_accessor( ElOp<spacedim>(this->dim_, {1}, this->n_columns_, false, nullptr, &side_reinit::ptop_weights) );
    this->n_columns_ += weights_side.n_comp();
    ElOp<spacedim> &JxW_side = this->add_accessor( ElOp<spacedim>(this->dim_, {1}, this->n_columns_, false, nullptr, &side_reinit::ptop_JxW, &weights_side) );
    this->n_columns_ += JxW_side.n_comp();
}

} // closing namespace FeSide

template class PatchPointValues<3>;
template class ElOp<3>;
template class FeBulk::PatchPointValues<3>;
template class FeSide::PatchPointValues<3>;
