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
 * Implementation of EvalSubset methods.
 */

EvalSubset::EvalSubset(std::shared_ptr<EvalPoints> eval_points, unsigned int n_permutations, unsigned int points_per_side)
: eval_points_(eval_points), subset_index_(eval_points_->n_subsets()), n_permutations_(n_permutations) {
    if (n_permutations_==0) { //bulk subset
        n_sides_ = 0;
        perm_indices_ = nullptr;
    } else {
        n_sides_ = eval_points->point_dim()+1;
        perm_indices_ = new unsigned int** [n_sides_];
        for (unsigned int i_side=0; i_side<n_sides_; ++i_side) {
            perm_indices_[i_side] = new unsigned int* [n_permutations_];
            for (unsigned int i_perm=0; i_perm<n_permutations_; ++i_perm) {
                perm_indices_[i_side][i_perm] = new unsigned int [points_per_side];
            }
        }
    }
}

EvalSubset::~EvalSubset() {
    /*if (perm_indices_!=nullptr) {
        for (unsigned int i_side=0; i_side<n_sides_; ++i_side) {
            for (unsigned int i_perm=0; i_perm<n_permutations_; ++i_perm) {
                delete perm_indices_[i_side][i_perm];
            }
            delete perm_indices_[i_side];
        }
        delete perm_indices_;
    }*/
}

Range< BulkPoint > EvalSubset::points(const DHCellAccessor &cell) const {
    ASSERT_EQ(n_sides_, 0).error("Method points with DHCellAccessor argument must be call for bulk subset!\n");
    if (cell.element_cache_index() == ElementCacheMap::undef_elem_idx)
        THROW( ExcElementNotInCache() << EI_ElementIdx(cell.elm_idx()) );

    auto bgn_it = make_iter<BulkPoint>( BulkPoint(cell, *this, eval_points_->subset_begin(subset_index_)) );
    auto end_it = make_iter<BulkPoint>( BulkPoint(cell, *this, eval_points_->subset_end(subset_index_)) );
    return Range<BulkPoint>(bgn_it, end_it);
}

Range< SidePoint > EvalSubset::points(const DHCellSide &cell_side) const {
    ASSERT_GT(n_sides_, 0).error("Method points with DHCellSide argument must be call for side subset!\n");
    if (cell_side.cell().element_cache_index() == ElementCacheMap::undef_elem_idx)
        THROW( ExcElementNotInCache() << EI_ElementIdx(cell_side.cell().elm_idx()) );

    unsigned int begin_idx = eval_points_->subset_begin(subset_index_);
    unsigned int end_idx = eval_points_->subset_end(subset_index_);
    unsigned int points_per_side = (end_idx - begin_idx) / this->n_sides();
    auto bgn_it = make_iter<SidePoint>( SidePoint(cell_side, *this, 0 ) );
    auto end_it = make_iter<SidePoint>( SidePoint(cell_side, *this, points_per_side ) );
    return Range<SidePoint>(bgn_it, end_it);
}



/******************************************************************************
 * Implementation of BulkPoint methods
 */


/******************************************************************************
 * Implementation of SidePoint methods
 */
