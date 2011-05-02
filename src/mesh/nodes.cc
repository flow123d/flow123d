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
 * @ingroup mesh
 * @brief  Class nodes
 *
 */

#if 0
#include "mesh/nodes.hh"

Node::Node() {
    coordinates = new double[3];

    n_elements = NDEF;

    scalar = 0.0;
    //conc = NULL;
    aux = NDEF;
    faux = 0.0;
}

void Node::set(double x, double y, double z) {
    coordinates[0] = x;
    coordinates[1] = y;
    coordinates[2] = z;
}

double Node::getX() {
    return coordinates[0];
}

double Node::getY() {
    return coordinates[1];
}

double Node::getZ() {
    return coordinates[2];
}

/**
 * Calculation of distance between nodes (given node and this node).
 *
 */
double Node::distance(Node* node) {
    ASSERT(!(node == NULL), "NULL as argument of method distance(TNode*)\n");

    double distance;

    distance = sqrt(
            (node->getX() - getX()) * (node->getX() - getX())
            + (node->getY() - getY()) * (node->getY() - getY())
            + (node->getZ() - getZ()) * (node->getZ() - getZ())
            );

    return distance;
}
#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
