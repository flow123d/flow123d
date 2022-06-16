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

EdgeIntegral::EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim, uint i_subset)
: BaseIntegral(eval_points, dim),
  subset_index_(i_subset)
{

    begin_idx_ = eval_points_->subset_begin(dim_, subset_index_);
    uint end_idx = eval_points_->subset_end(dim_, subset_index_);
    n_sides_ = dim + 1;
    //DebugOut() << "begin: " << begin_idx_ << "end: " << end_idx;
    n_points_per_side_ = (end_idx - begin_idx_) / n_sides();
    //DebugOut() << "points per side: " << n_points_per_side_;

}

EdgeIntegral::~EdgeIntegral() {
}


/******************************************************************************
 * Implementation of CouplingIntegral methods
 */

CouplingIntegral::CouplingIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral)
 : BaseIntegral(edge_integral->eval_points(), edge_integral->dim()),
   edge_integral_(edge_integral), bulk_integral_(bulk_integral)
{
    ASSERT_EQ(edge_integral->dim()-1, bulk_integral->dim());
}

CouplingIntegral::~CouplingIntegral()
{
    edge_integral_.reset();
    bulk_integral_.reset();
}



/******************************************************************************
 * Implementation of BoundaryIntegral methods
 */

BoundaryIntegral::BoundaryIntegral(std::shared_ptr<EdgeIntegral> edge_integral, std::shared_ptr<BulkIntegral> bulk_integral)
 : BaseIntegral(edge_integral->eval_points(), edge_integral->dim()),
   edge_integral_(edge_integral), bulk_integral_(bulk_integral)
{}

BoundaryIntegral::~BoundaryIntegral()
{
    edge_integral_.reset();
}
