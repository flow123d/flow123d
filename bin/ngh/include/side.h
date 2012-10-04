#ifndef sideH
#define sideH

#include "node.h"
#include "edge.h"
#include "shape.h"
#include "element.h"

#define FOR_SIDES(pSide)   for(;(pSide)!=NULL;(pSide)=(pSide)->getNext())
#define FOR_SIDE_NODES(n)  for(int iiSide=0; iiSide<n; iiSide++)

class TSide {
private:
    static int numberInstance;

    int id;
    int lnum; // Local # of side in element

    int n_nodes;
    TNode** node;

    TElement* ele;

    bool aux;

    TShape shape;

    TSide* prev;
    TSide* next;

    TEdge* edg;

    int generateId();
    void Initialize();

public:
    TSide(TElement*, int, TNode*);
    TSide(TElement*, int, TNode*, TNode*);
    TSide(TElement*, int, TNode*, TNode*, TNode*);
    ~TSide();

    int GetId();
    int GetLnum();
    int GetNNodes();
    bool AreBBNeighbours(TSide*);
    bool AreVBNeighbours(TElement*);

    void setAux(bool);
    bool getAux();

    void setPrev(TSide*);
    void setNext(TSide*);

    TSide* getNext();

    void setEdg(TEdge*);
    TEdge* getEdg();

    TElement* getElement();
    TNode* getNode(int);

    static int getNumInstances() {
        return TSide::numberInstance;
    }
};

#endif
