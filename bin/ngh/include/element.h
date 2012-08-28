#ifndef elementH
#define elementH

#include "shape.h"

#include "node.h"
#include "side.h"

#include "abscissa.h"
#include "triangle.h"
#include "tetrahedron.h"

#define FOR_ELEMENT_NODES(i,j)           for((j)=0;(j)<(i)->GetNNodes();(j)++)
#define FOR_ELEMENT_SIDES(i,j)           for((j)=0;(j)<(i)->GetNSides();(j)++)
#define FOR_ELEMENTS(ele, firstElement)  for((ele)=firstElement;(ele)!=NULL;(ele)=(ele)->getNext())

class TMesh;
class TProblem;

class TElement {
private:
    static int numberInstance;
    int id;

    int label;

    TShape type;

    int dim;

    TElement* prev;
    TElement* next;

    int n_nodes;
    TNode** node;

    int n_sides;
    TSide** side;

    int n_tags;

    int mid;
    int rid;

    double metrics;

    bool metrics_is_actual;

    TAbscissa* abscissa;
    TTriangle* triangle;
    TTetrahedron* tetrahedron;

    int generateId();

    void Allocate();
    void Initialize();

public:
    TElement(int, int, int, int*, TNode**);
    ~TElement();

    int GetNNodes();
    TNode* getNode(int);

    int GetNSides();

    int getLabel();
    int getNodeLabel(int);

    TShape GetType();

    void CreateSides(TMesh*);
    void CreateGeometry();

    void Test();

    void ComputeMetrics();
    double GetMetrics();

    bool AreBBNeighbours(TElement*);
    bool AreVBNeighbours(TElement*);
    bool AreVVNeighbours(TProblem*, TElement*, double&);

    void setSide(int, TSide*);

    void setPrev(TElement*);
    void setNext(TElement*);
    TElement* getNext();

    static int getNumInstances() {
        return TElement::numberInstance;
    }
};

#endif

