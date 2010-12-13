#ifndef PPFCS_H
#define PPFCS_H

struct Transport;

#define FOR_ELEMENTCUT(i)        for((i)=transport->fsec->elc;(i)!=NULL;(i)=(i)->next)

//=============================================================================
// STRUCTURE OF THE FLOW-CROSS-SECTION
//=============================================================================
struct FSection
{
		char *fcs_params;
		double      eqn[4] ;    // coef. of x,y,z,1
        int         axis_output;
        int         n_elm;    //  number of FCS elements
        struct  ElementCut    *elc;
};
//=============================================================================
// STRUCTURE OF THE ELEMENT CUT
//=============================================================================
struct ElementCut
{
        ElementIter element;
        int             type; // 1 - side break, 2 - element break
        int             sid;
        int             n_point;
        double          *point[4];     // coord.
        double          cutflux;
        struct ElementCut *next;
        struct ElementCut *prev;
};

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
