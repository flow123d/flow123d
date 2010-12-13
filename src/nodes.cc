/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Class nodes
 *
 */

#include "nodes.h"

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
//-----------------------------------------------------------------------------
// vim: set cindent:
