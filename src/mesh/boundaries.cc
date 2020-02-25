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
 * @file    boundaries.cc
 * @ingroup mesh
 * @brief   Boundary conditions
 */

// #include "system/system.hh"
#include "mesh/boundaries.h"
#include "mesh/mesh_data.hh"


flow::VectorId<unsigned int> Boundary::id_to_bcd;


Boundary::Boundary()
: boundary_data_(nullptr)
{}

Boundary::Boundary(BoundaryData* boundary_data)
: boundary_data_(boundary_data)
{}

//-----------------------------------------------------------------------------
// vim: set cindent:

