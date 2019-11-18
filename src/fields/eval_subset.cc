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

const EvalSubset EvalSubset::dummy_subset = EvalSubset();

EvalSubset::EvalSubset(std::shared_ptr<EvalPoints> eval_points)
: eval_points_(eval_points), block_indices_(1), n_sides_(0) {
	block_indices_[0] = eval_points_->n_blocks() - 1;
}

EvalSubset::EvalSubset(std::shared_ptr<EvalPoints> eval_points, unsigned int n_permutations)
: eval_points_(eval_points), block_indices_(n_permutations), n_sides_(eval_points->point_dim()+1) {
    for (unsigned int i=0, block_idx=eval_points_->n_blocks() - 1; i<n_permutations; ++i, ++block_idx)
    	block_indices_[i] = block_idx;
}

Range< BulkPoint > EvalSubset::points(const DHCellAccessor &cell) const {
    ASSERT_EQ(n_sides_, 0).error("Method points with DHCellAccessor argument must be call for bulk subset!\n");
    if (cell.element_cache_index() == ElementCacheMap::undef_elem_idx)
        THROW( ExcElementNotInCache() << EI_ElementIdx(cell.elm_idx()) );

    auto bgn_it = make_iter<BulkPoint>( BulkPoint(cell, *this, eval_points_->block_idx(block_indices_[0])) );
    auto end_it = make_iter<BulkPoint>( BulkPoint(cell, *this, eval_points_->block_idx(block_indices_[0]+1)) );
    return Range<BulkPoint>(bgn_it, end_it);
}

Range< SidePoint > EvalSubset::points(const DHCellSide &cell_side) const {
    ASSERT_GT(n_sides_, 0).error("Method points with DHCellSide argument must be call for side subset!\n");
    if (cell_side.cell().element_cache_index() == ElementCacheMap::undef_elem_idx)
        THROW( ExcElementNotInCache() << EI_ElementIdx(cell_side.cell().elm_idx()) );

    unsigned int perm_idx = cell_side.element()->permutation_idx( cell_side.side_idx() );
    unsigned int begin_idx = eval_points_->block_idx(block_indices_[perm_idx]);
    unsigned int end_idx = eval_points_->block_idx(block_indices_[perm_idx+1]);
    unsigned int points_per_side = (end_idx - begin_idx) / this->n_sides();
    auto bgn_it = make_iter<SidePoint>( SidePoint(cell_side, *this, begin_idx + cell_side.side_idx() * points_per_side ) );
    auto end_it = make_iter<SidePoint>( SidePoint(cell_side, *this, begin_idx + (cell_side.side_idx()+1) * points_per_side ) );
    return Range<SidePoint>(bgn_it, end_it);
}



/******************************************************************************
 * Implementation of BulkPoint methods
 */


/******************************************************************************
 * Implementation of SidePoint methods
 */
SidePoint SidePoint::permute(DHCellSide edg_side) const {
	unsigned int begin_idx = subset_.eval_points()->block_idx(subset_.get_block_idx(permutation_idx_));
	unsigned int points_per_side = (subset_.eval_points()->block_idx(subset_.get_block_idx(permutation_idx_+1)) - begin_idx) / subset_.n_sides();
	unsigned int point_on_side = local_point_idx_ - (begin_idx + cell_side_.side_idx() * points_per_side);

    unsigned int edg_perm_idx = edg_side.element()->permutation_idx( edg_side.side_idx() );
    unsigned int edg_begin_idx = subset_.eval_points()->block_idx( subset_.get_block_idx(edg_perm_idx) );
	unsigned int edg_loc_point_idx = edg_begin_idx + edg_side.side_idx() * points_per_side + point_on_side;
    return SidePoint(edg_side, subset_, edg_loc_point_idx);
}
