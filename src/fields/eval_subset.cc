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


/******************************************************************************
 * Implementation of EvalSubset methods.
 */

const EvalSubset EvalSubset::dummy_subset = EvalSubset();

const std::vector<int> &EvalSubset::get_point_indices(unsigned int permutation) const {
    ASSERT_LT( permutation, point_indices_.size() );
    return point_indices_[permutation];
}

Range< BulkPoint > EvalSubset::points(const DHCellAccessor &cell) const {
	ASSERT_EQ(n_sides_, 0).error("Method points with DHCellAccessor argument must be call for bulk subset!\n");
    auto bgn_it = make_iter<BulkPoint>( BulkPoint(cell, *this, 0) );
    auto end_it = make_iter<BulkPoint>( BulkPoint(cell, *this, point_indices_[0].size()) );
    return Range<BulkPoint>(bgn_it, end_it);
}

Range< SidePoint > EvalSubset::points(const DHCellSide &cell_side) const {
	ASSERT_GT(n_sides_, 0).error("Method points with DHCellSide argument must be call for side subset!\n");
	unsigned int perm_idx = cell_side.element()->permutation_idx( cell_side.side_idx() );
    unsigned int points_per_side = this->point_indices_[perm_idx].size() / this->n_sides();
    auto bgn_it = make_iter<SidePoint>( SidePoint(cell_side, *this, cell_side.side_idx()*points_per_side ) );
    auto end_it = make_iter<SidePoint>( SidePoint(cell_side, *this, (cell_side.side_idx()+1)*points_per_side ) );
    return Range<SidePoint>(bgn_it, end_it);
}



/******************************************************************************
 * Implementation of BulkPoint methods
 */
arma::vec BulkPoint::loc_coords() const {
    return this->eval_points().local_point( this->point_set_idx() );
}



/******************************************************************************
 * Implementation of SidePoint methods
 */
arma::vec SidePoint::loc_coords() const {
    return this->eval_points().local_point( this->point_set_idx() );
}
