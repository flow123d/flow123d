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
 * @file    neighbours.cc
 * @ingroup mesh
 * @brief   Initialize neighbouring
 */

#include "system/system.hh"
#include "neighbours.h"
#include "mesh/mesh.h"

//=============================================================================
// READ DATA OF ALL NEIGHBOURS
//=============================================================================
void read_neighbour_list(Mesh* mesh)
{

}

Neighbour::Neighbour()
: edge_idx_(-1)
{}

void Neighbour::reinit(ElementIter ele, unsigned int edg_idx)
{
    element_=ele;
    edge_idx_=edg_idx;
}


//-----------------------------------------------------------------------------
// vim: set cindent:
