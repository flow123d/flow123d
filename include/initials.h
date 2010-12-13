#ifndef INITIALS_H
#define INITIALS_H

class Mesh;

//=============================================================================
// STRUCTURE OF THE INITIAL CONDITION
//=============================================================================
struct Initial
{
    // Data readed from initial conditions files
    int      id;        // Id number of condition
    int      eid;           // Id number of element where prescribed
    double   epress;    // Press in the element's center
    int      n_sides;       // # of sides of the element readed from file
    double   *spress;       // Presses on the sides of the element
    double   *sflux;    // Fluxes via sides of the element
    int      n_neighs_vv;   // # of neighbours, V-V type (noncomp.)
    int  *hid;      // ID numbers of neighbour
    // Topology of the mesh
    // List
    struct Initial *prev;   // Previous initial in the list
    struct Initial *next;   // Next initial in the list
    // Misc
    int  aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
};

#define FOR_INITIALS(i)   for((i)=mesh->initial;(i)!=NULL;(i)=(i)->next)

void read_initial_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
