/*!
 *
 * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    mesh_data.hh
 * @brief   Internal mesh data classes.
 */

#ifndef MESH_DATA_H
#define MESH_DATA_H



class SideIter;

class EdgeData{
    public:
        EdgeData()
        : n_sides(0), side_(nullptr)
        {};
        // Topology of the mesh
        unsigned int  n_sides;   // # of sides of edge
        SideIter *side_; // sides of edge (could be more then two e.g. 1D mesh in 2d space with crossing )
};




class Mesh;

class BoundaryData
{
public:
    static const unsigned int undef_idx=-1;

    BoundaryData()
    : edge_idx_(undef_idx),
      bc_ele_idx_(undef_idx),
      mesh_(nullptr)
    {};

    // Topology of the mesh
    unsigned int    edge_idx_;    // more then one side can be at one boundary element
    unsigned int    bc_ele_idx_;  // in near future this should replace Boundary itself, when we remove BC data members
    Mesh *mesh_;

};

#endif