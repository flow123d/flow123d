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

#ifndef NODE_H
#define NODE_H

#include <mesh_types.hh>
#include <point.h>

using namespace dealii;

/**
 * Class of node.
 * First approach in turning to class.
 */
class Node {
private:
    /** Coordinates of node */
    double* coordinates;
    Point<3>    point;

public:
    Node();

    void set(double, double, double);

    double getX();
    double getY();
    double getZ();

    double distance(Node*);

    int id;

    //--------------------------------------------------------------------------
    // Old data - adepts to reduce ...
    //--------------------------------------------------------------------------

    // Topology
    int n_elements; // # of elms connected by the node ( used only in transport )
    ElementIter *element; // List of elements

    // following  is used only by interpolation function
    // postprocess.c void make_node_scalar(Mesh* mesh)
    // which should be rewrittento be able interpolate arbitrary data
    // Results
    double scalar; // Scalar quantity (pressure/piez. head)

    // Misc
    int aux; // Auxiliary flag
    double faux; // Auxiliary number
};

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
