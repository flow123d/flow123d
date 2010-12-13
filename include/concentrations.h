#ifndef CONCENTRATIONS_H
#define CONCENTRATIONS_H

#include "mesh.h"

//=============================================================================
// STRUCTURE OF THE CONCENTRATION
//=============================================================================
struct Concentration
{
    // Data readed from concentration files
    int id;        // Id number of condition
    int eid;       // ID number of element where prescribed
    int n_subst;   // # of substances
    double *conc;         // Values of concentrations
    // List
    struct Concentration *prev; // Previous concentration in the list
    struct Concentration *next; // Next concentration in the list
    // Misc
    int aux;       // Auxiliary flag
    double faux;      // Auxiliary number
};

#define FOR_CONCENTRATIONS(i)   for((i)=mesh->concentration;(i)!=NULL;(i)=(i)->next)

void read_concentration_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
