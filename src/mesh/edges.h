/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief ???
 *
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
