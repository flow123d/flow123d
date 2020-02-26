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
 * @file    boundaries.h
 * @brief   
 */

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "mesh/region.hh"    // for Region
#include "mesh/sides.h"      // for SideIter
#include "mesh/mesh.h"       // for Mesh
#include "mesh/accessors.hh"  // for ElementAccessor
#include "mesh/mesh_data.hh"  // for BoundaryData
#include "system/system.hh"  // for MessageType::Err, xprintf

class Element;


/**
 * Setting boundary conditions should have two staps.
 * 1) Denote by numbers segments of mesh boundary. Possibly every side can be boundary.
 * 2) Assign particular type and values of BC on every boundary segment.
 *
 * So in future Boundary should keep only side and segment and there should be
 * one Boundary for every external side. Side is external either when it does not
 * neighbor with another element or when it belongs to an segment.
 */

//=============================================================================
// STRUCTURE OF THE BOUNDARY CONDITION
//=============================================================================
class Boundary
{
public:
    Boundary();
    Boundary(BoundaryData* boundary_data);

    inline bool is_valid() const {
        return boundary_data_ != nullptr;
    }

    inline Edge edge()
    {
        ASSERT_DBG(is_valid());
        return boundary_data_->mesh_->edge(boundary_data_->edge_idx_);
    }

    inline Element * element()
    {
        ASSERT_DBG(is_valid());
        return &( boundary_data_->mesh_->element_vec_[boundary_data_->bc_ele_idx_] );
    }

    inline Region region()
    { return element_accessor().region(); }

    inline ElementAccessor<3> element_accessor()
    {
        ASSERT_DBG(is_valid());
        return boundary_data_->mesh_->element_accessor(boundary_data_->bc_ele_idx_);
    }

    inline Mesh* mesh()
    {
        ASSERT_DBG(is_valid());
        return boundary_data_->mesh_;
    }

    inline uint edge_idx()
    {
        ASSERT_DBG(is_valid());
        return boundary_data_->edge_idx_;
    }

    inline uint bc_ele_idx()
    {
        ASSERT_DBG(is_valid());
        return boundary_data_->bc_ele_idx_;
    }

private:
    BoundaryData* boundary_data_;
};
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
