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
 * @file    edges.h
 * @brief   
 */

#ifndef MAKE_EDGES_H
#define MAKE_EDGES_H

#include "mesh/mesh.h"

//=============================================================================
// STRUCTURE OF THE EDGE OF THE MESH
//=============================================================================
class Edge
{
public:
    /// Minimalistic default constructor.
    Edge();
    inline SideIter side(const unsigned int i) const {
        return side_[i];
    }

    // Topology of the mesh
    int  n_sides;   // # of sides of edge
    SideIter *side_; // sides of edge (could be more then two e.g. 1D mesh in 2d space with crossing )

};

#define FOR_EDGE_SIDES(i,j) for((j)=0;(j)<(i)->n_sides;(j)++)


#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
