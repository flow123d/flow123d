#ifndef meshH
#define meshH

#include <map>

#include "element.h"
#include "node.h"
#include "side.h"
#include "edge.h"
#include "neighbour.h"

class TMesh {
private:
    int n_sides;
    int n_edges;
    int n_neighbours;

    TNode* firstNode;
    TNode* lastNode;

    TElement* firstElement;
    TElement* lastElement;

    TSide* firstSide;
    TSide* lastSide;

    TEdge* firstEdge;
    TEdge* lastEdge;

    TNeighbour* first_ngh;
    TNeighbour* last_ngh;

    std::map<int, TNode*> nodeLabelMap;
    std::map<int, TElement*> elementLabelMap;

public:
    TMesh();
    ~TMesh();

    void AddNode(TNode*);
    TNode* getNodeLabel(int);

    void AddElement(TElement*);
    TElement* getElementLabel(int);
    int getNumElements();

    void AddSide(TSide*);
    void AddEdge(TEdge*);
    void AddNeighbour(TNeighbour*);

    void NodeToElement();
    void SideToNode();

    void CreateSides();
    void CreateEdges();

    int GetNNeighs();

    void CreateGeometry();
    void TestElements();

    TElement* getFirstElement();
    TSide* getFirstSide();
    TNeighbour* getFirstNGH();

    void CreateNeighboursBB();
    void CreateNeighboursVB();
};
#endif
