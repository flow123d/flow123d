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

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

//#include "mesh/mesh.h"
#include "mesh/sides.h"
#include "mesh/edges.h"
#include "system/sys_vector.hh"



/**
 * Setting boundary conditions should have two staps.
 * 1) Denote by numbers segments of mesh boundary. Possibly every side can be boundary.
 * 2) Assign particular type and values of BC on every boundary segment.
 *
 * So in future Boundary should keep only side and segment and there should be
 * one Boundary for every external side. Side is external either when it does not
 * neighbor with another element or when it belongs to an segment.
 */

class Element;

//=============================================================================
// STRUCTURE OF THE BOUNDARY CONDITION
//=============================================================================
class Boundary
{
public:
    /**
     * temporary solution for old type BCD.
     * Transport BCD refers through IDs to flow BCD, so we have to
     * store positions of Flow BCD items somewhere.
     */
    static flow::VectorId<unsigned int> id_to_bcd;

    Boundary();

    /**
     * Can not make this inline now.
     */
    Edge * edge();

    Element * element();

    Region region() {
        return element()->region();
    }

    ElementAccessor<3> element_accessor();


    inline SideIter side() {
        if (edge()->n_sides != 1) xprintf(Err, "Using side method for boundary, but there is boundary with multiple sides.\n");
        return edge()->side_[0];
    }

    // Topology of the mesh
    unsigned int    edge_idx_;    // more then one side can be at one boundary element
    unsigned int    bc_ele_idx_;  // in near future this should replace Boundary itself, when we remove BC data members
    Mesh *mesh_;

};
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
