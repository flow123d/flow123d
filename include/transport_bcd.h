#ifndef TRANSPORT_BCD_H
#define TRANSPORT_BCD_H

class Mesh;

//=============================================================================
// STRUCTURE OF THE TRANSPORT_BCD
//=============================================================================
struct Transport_bcd {
    // Data readed from transport bounary conditions files
    int id; // Id number of condition
    int bid; // ID number of boundary condition where prescribed
    int n_subst; // # of substances
    double *conc; // Values of concentrations

    // List
    struct Transport_bcd *prev; // Previous transport boundary in the list
    struct Transport_bcd *next; // Next transport boundary in the list

    // Misc
    int aux; // Auxiliary flag
    double faux; // Auxiliary number
};

#define FOR_TRANSPORT_BCDS(i)   for((i)=mesh->transport_bcd;(i)!=NULL;(i)=(i)->next)

void read_transport_bcd_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
