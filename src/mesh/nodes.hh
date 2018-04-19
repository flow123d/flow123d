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
 * @file    nodes.hh
 * @brief   Nodes of a mesh.
 */

#ifndef NODE_H
#define NODE_H

#include "system/global_defs.h"
#include "mesh/mesh_types.hh"
#include <armadillo>
#include "system/armor.hh"




/**
 * Class of node.
 * First approach in turning to class.
 */
class Node {
private:
    /// Node point in 3D space.
    armor::vec<3> coordinates;

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
        {coordinates(0)=x; coordinates(1)=y; coordinates(2)=z;}

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
    /*
    inline arma::vec3 &point()
    { return coordinates.arma(); }
    */
    
    inline void set_point(double x, double y, double z) {
        coordinates(0) = x;
        coordinates(1) = y;
        coordinates(2) = z;
    }
    
    inline void set_point(arma::vec3 a) {
        coordinates = a;
    }
    
    inline arma::vec3 point() const
        { return coordinates.arma(); }
        
    inline const armor::vec<3> & point_armor() const {
        return coordinates;
    }
    
    inline armor::vec<3> & point_armor() {
        return coordinates;
    }

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


    // Misc
    int aux; // Auxiliary flag
};

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
