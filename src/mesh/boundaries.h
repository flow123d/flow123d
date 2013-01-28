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

#include "mesh/mesh.h"
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


//=============================================================================
// STRUCTURE OF THE BOUNDARY CONDITION
//=============================================================================
class Boundary
{
public:
    Boundary()
    : group(0), type(2), flux(0.0)
    {}
    /**
     * temporary solution for old type BCD.
     * Transport BCD refers through IDs to flow BCD, so we have to
     * store positions of Flow BCD items somewhere.
     */
    static flow::VectorId<unsigned int> id_to_bcd;

    inline ElementIter get_bc_element_iter() {
        return bc_element_;
    }

    // Data readed from boundary conditions files (REMOVE)
    int      type;      // Type of boundary condition
    double   scalar;    // Scalar - for Dirichlet's or Newton's type
    double   flux;      // Flux - for Neumann's type or source
    double   sigma;     // Sigma koef. - for Newton's type

    int      group;     // Group of condition
    // Topology of the mesh
    SideIter side;      // side, where prescribed
    ElementIter bc_element_;  // in near future this should replace Boundary itself, when we remove BC data members

};
#define DIRICHLET   1
#define NEUMANN     2
#define NEWTON      3


void read_boundary(Mesh*, const string &boundary_filename);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
