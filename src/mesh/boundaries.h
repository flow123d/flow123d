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


class Element;
class Mesh;
class Edge;
class Region;
class SideIter;
namespace flow { template <class T> class VectorId; }
template <int spacedim> class ElementAccessor;



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

    Region region();

    ElementAccessor<3> element_accessor();

    SideIter side();

    // Topology of the mesh
    unsigned int    edge_idx_;    // more then one side can be at one boundary element
    unsigned int    bc_ele_idx_;  // in near future this should replace Boundary itself, when we remove BC data members
    Mesh *mesh_;
};
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
