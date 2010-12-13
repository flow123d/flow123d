#ifndef MAKE_NEIGHBOURS_H
#define MAKE_NEIGHBOURS_H

#include <mesh_types.hh>

struct Problem;
class Mesh;
struct Element;
struct Neighbour;
struct Side;
struct Edge;

//=============================================================================
// STRUCTURE OF THE NEIGHBOUR
//=============================================================================
typedef struct Neighbour
{
    int   id;           // Global id number
    int   type;         // Type
    int   n_sides;      // # of neighbouring sides
    int   n_elements;   // # of neighbouring elements
    //int   n_coefs;      // # of coefficients
    int  *sid;      // id numbers of neighbouring sides
    int  *eid;      // id numbers of neighbouring elements
    double flux;           // flux for VV neigh. from e1 into e2 is negative
                                // from e2 into e1 is positive
    //double *coef;       // coefficients of neighbouring
    double  sigma;      // transition coefficient
    double geom_factor; // fraction of the lower dimension element for VV ngh
    struct Edge *edge;  // edge (set of neighbouring sides)
    struct Side **side; // neighbouring sides
    ElementIter *element;  // neighbouring elements
                               // for VB  - element[0] is element of lower dimension
    char *line;     // Type specific part of line of the input file
    // List
    struct Neighbour *prev; // Previous neighbour in the list
    struct Neighbour *next; // Next neighbour in the list
    // Misc
    int  aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
} Neighbour;

// Input neigbouring codes
#define BB_E         10     // two elements of same dim specified by eid
#define BB_EL        11     // two elements ... specified by explicit eid and sid
#define VB_ES        20     // compatible
#define VV_2E        30     // noncompatible

#define FOR_NEIGHBOURS(i)   for((i)=mesh->neighbour;(i)!=NULL;(i)=(i)->next)
#define FOR_NEIGH_ELEMENTS(i,j) for((j)=0;(j)<(i)->n_elements;(j)++)
#define FOR_NEIGH_SIDES(i,j)    for((j)=0;(j)<(i)->n_sides;(j)++)

void read_neighbour_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
