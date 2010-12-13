#ifndef SOURCES_H
#define SOURCES_H

#include <mesh_types.hh>

class Mesh;

//=============================================================================
// STRUCTURE OF THE SOURCE
//=============================================================================

struct Source {
    int id;
    int type;
    int eid;
    double density;
    ElementIter element;
    struct Source* prev;
    struct Source* next;
    int aux;
    double faux;
};


#define FOR_SOURCES(i)      for((i)=mesh->source;(i)!=NULL;(i)=(i)->next)

void read_source_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:

