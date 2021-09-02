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

EdgeIntegral::EdgeIntegral(std::shared_ptr<EvalPoints> eval_points, unsigned int dim)
:  BaseIntegral(eval_points, dim), subset_index_(eval_points_->n_subsets(dim)) {
    n_sides_ = dim_+1;
}

EdgeIntegral::~EdgeIntegral() {
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
