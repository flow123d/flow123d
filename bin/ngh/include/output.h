#ifndef outputH
#define outputH

#include "config.h"
#include "mesh.h"
#include "neighbour.h"

void Output(char*, TMesh*);
void WriteNeighboursBB(TNeighbour*, FILE*);
void WriteNeighboursVB(TNeighbour*, FILE*);
void WriteNeighboursVV(TNeighbour*, FILE*);

#endif

