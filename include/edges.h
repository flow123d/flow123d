#ifndef MAKE_EDGES_H
#define MAKE_EDGES_H

#include "mesh.h"

//=============================================================================
// STRUCTURE OF THE EDGE OF THE MESH
//=============================================================================
typedef struct Edge
{
    // Basic
    int  id;        // Id # of the edge
    // Topology of the mesh
    int  n_sides;   // # of sides of edge
    struct Side **side; // sides of edge (could be more then two e.g. 1D mesh in 2d space with crossing )
    struct Neighbour *neigh_vb; // "Compatible" neighbouring
    struct Neighbour *neigh_bb; // ??? this is what
    // List
    struct Edge *prev;  // Previous edge in the list
    struct Edge *next;  // Next edge in the list
    // Matrix
    int  c_row;     // # of row in block C (and E and F) (MH)
    double  f_val;      // diagonal value  in block F
    double  f_rhs;      // rhs value
    // Misc
    int      aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
} Edge;

#define FOR_EDGES(i)        for((i)=mesh->edge;(i)!=NULL;(i)=(i)->next)
#define FOR_EDGE_SIDES(i,j) for((j)=0;(j)<(i)->n_sides;(j)++)

void make_edge_list(Mesh*);
void edge_calculation_mh(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
