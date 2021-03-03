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
 * @file    eval_subset.cc
 * @brief
 * @author  David Flanderka
 */

#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"


/******************************************************************************
 * Implementation of BaseIntegral methods
 */

BaseIntegral::~BaseIntegral()
{}


/******************************************************************************
 * Implementation of BulkIntegral methods
 */

BulkIntegral::~BulkIntegral()
{}


/******************************************************************************
 * Implementation of EdgeIntegral methods
 */

EdgeIntegral::EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, unsigned int n_permutations, unsigned int points_per_side)
:  BaseIntegral(eval_points, dim), subset_index_(eval_points_->n_subsets(dim)), n_permutations_(n_permutations) {
    n_sides_ = dim_+1;
    perm_indices_ = new unsigned int** [n_sides_];
    for (unsigned int i_side=0; i_side<n_sides_; ++i_side) {
        perm_indices_[i_side] = new unsigned int* [n_permutations_];
        for (unsigned int i_perm=0; i_perm<n_permutations_; ++i_perm) {
            perm_indices_[i_side][i_perm] = new unsigned int [points_per_side];
        }
    }
}

EdgeIntegral::~EdgeIntegral() {
    for (unsigned int i_side=0; i_side<n_sides_; ++i_side) {
        for (unsigned int i_perm=0; i_perm<n_permutations_; ++i_perm) {
            delete perm_indices_[i_side][i_perm];
        }
        delete perm_indices_[i_side];
    }
    delete perm_indices_;
}


/******************************************************************************
 * Implementation of CouplingIntegral methods
 */

CouplingIntegral::CouplingIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral)
 : BaseIntegral(edge_integral->eval_points(), edge_integral->dim()), edge_integral_(edge_integral), bulk_integral_(bulk_integral) {
    ASSERT_EQ_DBG(edge_integral->dim(), bulk_integral->dim());
}

CouplingIntegral::~CouplingIntegral() {
    edge_integral_.reset();
    bulk_integral_.reset();
}



/******************************************************************************
 * Implementation of BoundaryIntegral methods
 */

BoundaryIntegral::BoundaryIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral)
 : BaseIntegral(edge_integral->eval_points(), edge_integral->dim()), edge_integral_(edge_integral), bulk_integral_(bulk_integral) {}

BoundaryIntegral::~BoundaryIntegral() {
    edge_integral_.reset();
}


/******************************************************************************
 * Implementation of BulkPoint methods
 */


/******************************************************************************
 * Implementation of EdgePoint methods
 */

EdgePoint EdgePoint::point_on(DHCellSide edg_side) const {
    return EdgePoint(edg_side, elm_cache_map_, this->integral_, this->local_point_idx_);
}

unsigned int EdgePoint::eval_point_idx() const {
    return integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}


/******************************************************************************
 * Implementation of CouplingPoint methods
 */

BulkPoint CouplingPoint::lower_dim(DHCellAccessor cell_lower) const {
    return BulkPoint(cell_lower, elm_cache_map_, integral_->bulk_integral_.get(),
            this->eval_points()->subset_begin(cell_lower.dim(), integral_->get_subset_low_idx())+local_point_idx_);
}

unsigned int CouplingPoint::eval_point_idx() const {
    return integral_->edge_integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}


/******************************************************************************
 * Implementation of BoundaryPoint methods
 */

BulkPoint BoundaryPoint::point_bdr(ElementAccessor<3> bdr_elm) const {
    return BulkPoint(bdr_elm, elm_cache_map_, integral_->bulk_integral_.get(),
            this->eval_points()->subset_begin(bdr_elm.dim(), integral_->get_subset_low_idx())+local_point_idx_);
}

unsigned int BoundaryPoint::eval_point_idx() const {
    return integral_->edge_integral_->perm_idx_ptr(side_idx_, permutation_idx_, local_point_idx_);
}
