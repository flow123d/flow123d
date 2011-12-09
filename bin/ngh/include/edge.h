#ifndef edgeH
#define edgeH

#include <vector>

#define FOR_EDGES(i,mesh)    for((i)=mesh->firstEdge;(i)!=NULL;(i)=(i)->getNext())
#define FOR_EDGE_SIDES(i,j)  for((j)=0;(j)<(i)->GetNSides();(j)++)

class TSide;

class TEdge {
private:
    static int numberInstance;
    int id;

    TEdge* prev;
    TEdge* next;

    int n_sides;
    int* sid;
    TSide** side;

public:
    TEdge(std::vector<TSide*>);
    ~TEdge();

    int GetId();
    int GetNSides();

    int getSId(int);

    void setPrev(TEdge*);
    TEdge* getPrev();

    void setNext(TEdge*);
    TEdge* getNext();

    TSide* getSide(int);

    static int getNumInstances() {
        return TEdge::numberInstance;
    }
};

#endif
