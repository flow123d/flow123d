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
