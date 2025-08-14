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
 * @file    integral_points.cc
 * @brief
 * @author  David Flanderka
 */

#include "fem/integral_points.hh"
#include "fem/integral_acc.hh"



EdgePoint EdgePoint::point_on(const DHCellSide &edg_side) const {
    uint element_patch_idx = elm_cache_map_->position_in_cache(edg_side.element().idx());
    uint side_begin = integral_->side_begin(edg_side);
    return EdgePoint(BulkPoint(elm_cache_map_, element_patch_idx, local_point_idx_),
            integral_, side_begin);
}

//******************************************************************************
BulkPoint CouplingPoint::lower_dim(DHCellAccessor cell_lower) const {
    unsigned int i_elm = elm_cache_map_->position_in_cache(cell_lower.elm().idx());
    unsigned int i_ep = integral_->bulk_begin() + local_point_idx_;
    return BulkPoint(elm_cache_map_, i_elm, i_ep);
}



//******************************************************************************
BulkPoint BoundaryPoint::point_bdr(ElementAccessor<3> bdr_elm) const {
    unsigned int i_elm = elm_cache_map_->position_in_cache(bdr_elm.idx(), true);
    unsigned int i_ep = integral_->bulk_begin() + local_point_idx_;
    //DebugOut() << "begin:" << integral_->bulk_begin() << "iloc " << local_point_idx_;
    return BulkPoint(elm_cache_map_, i_elm, i_ep);
}


/******************************************************************************
 * Temporary implementations. Intermediate step in implementation of PatcFEValues.
 */
//
//unsigned int EdgePoint::side_idx() const {
//    return (this->side_begin_ - integral_->begin_idx_) / integral_->n_points_per_side_;
//}
//
//unsigned int CouplingPoint::side_idx() const {
//    return (this->side_begin_ - integral_->edge_integral_->begin_idx_) / integral_->edge_integral_->n_points_per_side_;
//}
//
//unsigned int BoundaryPoint::side_idx() const {
//    return (this->side_begin_ - integral_->edge_integral_->begin_idx_) / integral_->edge_integral_->n_points_per_side_;;
//}
