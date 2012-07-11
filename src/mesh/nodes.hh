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
 * @brief Nodes of a mesh.
 *
 */

#ifndef NODE_H
#define NODE_H

#include "system/global_defs.h"
#include "mesh/mesh_types.hh"
#include <armadillo>




/**
 * Class of node.
 * First approach in turning to class.
 */
class Node {
private:
    /// Node point in 3D space.
    arma::vec3 coordinates;

public:
    /**
     * Default constructor.
     */
    Node()
    //: element(NULL)
        {coordinates.zeros();}

    /**
     * Construct form given coordinates.
     *
     * Possibly there could be also constructor from a vector.
     */
    Node(double x, double y, double z)
    //: element(NULL)
        {coordinates(0)=x; coordinates(1)=y; coordinates(3)=z;}

    /**
     * Old getter methods. OBSOLETE.
     */
    inline double getX() const
        { return coordinates[0];}
    inline double getY() const
        {return coordinates[1];}
    inline double getZ() const
        {return coordinates[2];}


    /**
     * Getter method for nodal point. Can be used also for modification.
     */
    inline arma::vec3 &point()
    { return coordinates; }

    inline const arma::vec3 &point() const
        { return coordinates; }

    /**
     * Difference of two nodes is a vector.
     */
    inline arma::vec3 operator-(const Node &n2) const
    {
        return ( this->point() - n2.point() );
    }


    /**
     * Distance of two nodes.
     */
    inline double distance(const Node &n2) const
    {
        return norm(*this - n2, 2);
    }

    //--------------------------------------------------------------------------
    // Old data - adepts to remove ...
    //--------------------------------------------------------------------------

    //int id; // this is used in old output TODO: remove after application new output

    // Topology
    // int n_elements; // # of elms connected by the node ( used only in transport )
    // ElementIter *element; // List of elements

    // following  is used only by interpolation function
    // postprocess.c void make_node_scalar(Mesh* mesh)
    // which should be rewrittento be able interpolate arbitrary data
    // Results
    //double scalar; // Scalar quantity (pressure/piez. head)

    // Misc
    int aux; // Auxiliary flag
    //double faux; // Auxiliary number
};

/**
 * Binary operators (and functions) operating on nodes.
 */

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
