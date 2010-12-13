#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "mesh.h"

struct Boundary;
struct Transport_bcd;

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
    // Data readed from boundary conditions files (REMOVE)
    int      type;      // Type of boundary condition
    double   scalar;    // Scalar - for Dirichlet's or Newton's type
    double   flux;      // Flux - for Neumann's type or source
    double   sigma;     // Sigma koef. - for Newton's type

    int      group;     // Group of condition
    // Topology of the mesh
    struct Side *side;      // side, where prescribed
    struct Transport_bcd    *transport_bcd;  // transport boundary condition (REMOVE)

    // Misc
    int  aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
};
#define DIRICHLET   1
#define NEUMANN     2
#define NEWTON      3


void read_boundary(Mesh*);
void boundary_calculation_mh(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
