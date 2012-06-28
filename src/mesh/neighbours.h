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

#ifndef MAKE_NEIGHBOURS_H
#define MAKE_NEIGHBOURS_H

#include "mesh/mesh_types.hh"



struct Element;
struct Edge;
class SideIter;

//=============================================================================
// STRUCTURE OF THE NEIGHBOUR
//=============================================================================
class Neighbour
{
public:
    Neighbour();

    // side of the edge in higher dim. mesh
    inline SideIter side();

    // edge of lower dimensional mesh in VB neigh.
    inline Edge *edge();

    // element of higher dimension mesh in VB neigh.
    inline ElementIter element();

    int   type;         // Type
    int   n_sides;      // # of neighbouring sides
    int   n_elements;   // # of neighbouring elements

    int  *sid;      // id numbers of neighbouring sides
    int  *eid;      // id numbers of neighbouring elements
    double  sigma;      // transition coefficient
    struct Edge *edge_;  // edge (set of neighbouring sides)


    SideIter *side_; // neighbouring sides (THIS IS ONLY meaningfull member for VB )

    ElementIter *element_;  // neighbouring elements
                               // for VB  - element[0] is element of lower dimension
};

// Input neigbouring codes
#define BB_E         10     // two elements of same dim specified by eid
#define BB_EL        11     // two elements ... specified by explicit eid and sid
#define VB_ES        20     // compatible
#define VV_2E        30     // noncompatible



//void read_neighbour_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
