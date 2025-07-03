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

#include "fem/eval_subset.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"


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

EdgeIntegral::~EdgeIntegral() {
}


/******************************************************************************
 * Implementation of CouplingIntegral methods
 */

CouplingIntegral::~CouplingIntegral()
{
    edge_integral_.reset();
    bulk_integral_.reset();
}



/******************************************************************************
 * Implementation of BoundaryIntegral methods
 */

BoundaryIntegral::~BoundaryIntegral()
{
    edge_integral_.reset();
}


/******************************************************************************
 * Temporary implementations. Intermediate step in implementation of PatcFEValues.
 */

unsigned int EdgePoint::side_idx() const {
    return (this->side_begin_ - integral_->begin_idx_) / integral_->n_points_per_side_;
}

unsigned int CouplingPoint::side_idx() const {
    return (this->side_begin_ - integral_->edge_integral_->begin_idx_) / integral_->edge_integral_->n_points_per_side_;
}

unsigned int BoundaryPoint::side_idx() const {
    return (this->side_begin_ - integral_->edge_integral_->begin_idx_) / integral_->edge_integral_->n_points_per_side_;;
}
