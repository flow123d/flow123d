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
 * @file    integral_acc.cc
 * @brief
 * @author  David Flanderka
 */

#include "fem/integral_points.hh"
#include "fem/integral_acc.hh"
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
    internal_edge_.reset();
    internal_bulk_.reset();
}



/******************************************************************************
 * Implementation of BoundaryIntegral methods
 */

BoundaryIntegral::BoundaryIntegral(std::shared_ptr<EvalPoints> eval_points, Quadrature *quad, unsigned int dim)
 : BaseIntegral(eval_points, dim)
{
    switch (dim) {
    case 1:
        internal_bulk_ = eval_points->add_bulk_internal<0>(quad);
        internal_edge_ = eval_points->add_edge_internal<1>(quad);
        break;
    case 2:
        internal_bulk_ = eval_points->add_bulk_internal<1>(quad);
        internal_edge_ = eval_points->add_edge_internal<2>(quad);
        break;
    case 3:
        internal_bulk_ = eval_points->add_bulk_internal<2>(quad);
        internal_edge_ = eval_points->add_edge_internal<3>(quad);
        break;
    default:
        ASSERT(false).error("Should not happen!\n");
    }
}

BoundaryIntegral::~BoundaryIntegral()
{
    internal_edge_.reset();
    internal_bulk_.reset();
}
