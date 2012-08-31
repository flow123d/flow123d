#ifndef neighbourH
#define neighbourH

#include "edge.h"
#include "side.h"
#include "element.h"
#include "intersectionLocal.h"

#define FOR_NEIGHBOURS(i)        for((i)=mesh->getFirstNGH();(i)!=NULL;(i)=(i)->getNext())
#define FOR_NEIGH_ELEMENTS(i,j)  for((j)=0;(j)<(i)->GetNElements();(j)++)

class TEdge;
class TSide;
class TElement;
class TMesh;
class IntersectionLocal;

typedef enum NghTypes {
    nt_unknown = 0,
    nt_bb = 11,
    nt_vb = 20,
    nt_vv = 30
    //nt_vv = 31	//nekompatibilni + kompatibilni bb
    //nt_vv = 32	//nekompatibilni + kompatibilni vb
} TNghType;

class TNeighbour {
private:
    static int numberInstance;
    int id;

    TNghType type;

    //objekt typu IntersectionLocal
    IntersectionLocal* intersection;

    int n_sides;
    TSide** side;

    int n_elements;
    TElement** element;

    TNeighbour* prev;
    TNeighbour* next;

    TEdge* edge;

    int generateId();
    void Initialize();

public:
    TNeighbour(TEdge*);
    TNeighbour(TElement*, TSide*);
    TNeighbour(TElement*, TElement*, IntersectionLocal *);
    ~TNeighbour();

    int GetId();
    int GetNElements();
    //double GetCoef();
    TNghType GetType();

    TSide* getSide(int);
    TElement* getElement(int);
    IntersectionLocal *get_intersection() { return intersection; }

    void setPrev(TNeighbour*);
    TNeighbour* getPrev();
    void setNext(TNeighbour*);
    TNeighbour* getNext();

    static int getNumInstances() {
        return TNeighbour::numberInstance;
    }
};

#endif
