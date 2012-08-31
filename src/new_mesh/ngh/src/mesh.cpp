#include <iostream>
#include <vector>

#include "new_mesh/ngh/include/mesh.h"
#include "new_mesh/ngh/include/system.h"
#include "new_mesh/ngh/include/TElementFactory.h"

TMesh::TMesh() {
    firstNode = NULL;
    lastNode = NULL;

    firstElement = NULL;
    lastElement = NULL;

    firstSide = NULL;
    lastSide = NULL;

    firstEdge = NULL;
    lastEdge = NULL;

    first_ngh = NULL;
    last_ngh = NULL;

    n_sides = 0;
    n_edges = 0;
    n_neighbours = 0;
}

TMesh::~TMesh() {

    TElement* ele;

    FOR_ELEMENTS(ele, this->getFirstElement()) {
        //        delete ele;
    }
    elementLabelMap.clear();

    TNode* node;

    FOR_NODES(node, this) {
        //        delete node;
    }
    nodeLabelMap.clear();
}

void TMesh::AddNode(TNode* nod) {
    nodeLabelMap[nod->getLabel()] = nod;

    if (firstNode == NULL) {
        firstNode = nod;
        lastNode = nod;
    } else {
        lastNode->setNext(nod);
        nod->setPrev(lastNode);
        lastNode = nod;
    }

    return;
}

void TMesh::AddElement(TElement* ele) {
    elementLabelMap[ele->getLabel()] = ele;

    if (firstElement == NULL) {
        firstElement = ele;
        lastElement = ele;
    } else {
        lastElement->setNext(ele);
        ele->setPrev(lastElement);
        lastElement = ele;
    }

    return;
}

TNode* TMesh::getNodeLabel(int nodeLabel) {
    return nodeLabelMap[nodeLabel];
}

TElement* TMesh::getElementLabel(int elementLabel) {
    return elementLabelMap[elementLabel];
}

int TMesh::getNumElements() {
    return elementLabelMap.size();
}

TElement* TMesh::getFirstElement() {
    return firstElement;
}

TSide* TMesh::getFirstSide() {
    return firstSide;
}

TNeighbour* TMesh::getFirstNGH() {
    return first_ngh;
}

int TMesh::GetNNeighs() {
    return n_neighbours;
}

void TMesh::NodeToElement() {
    std::cout << "Assigning nodes to the elements and back... ";

    TElement* ele;

    FOR_ELEMENTS(ele, this->getFirstElement()) {
        int i;

        FOR_ELEMENT_NODES(ele, i) {
            TNode* node = ele->getNode(i);
            node->addElement(ele);
        }
    }
    std::cout << "OK.\n";
}

void TMesh::CreateSides() {
    std::cout << "Creating sides... ";

    TElement* ele;

    FOR_ELEMENTS(ele, this->getFirstElement()) {
        ele->CreateSides(this);
    }
    std::cout << n_sides << " sides created. OK.\n";
}

void TMesh::AddSide(TSide* sid) {
    n_sides++;

    if (firstSide == NULL) {
        firstSide = sid;
        lastSide = sid;
    } else {
        lastSide->setNext(sid);
        sid->setPrev(lastSide);
        lastSide = sid;
    }

    return;
}

void TMesh::AddEdge(TEdge* edg) {
    n_edges++;

    if (firstEdge == NULL) {
        firstEdge = edg;
        lastEdge = edg;
    } else {
        lastEdge->setNext(edg);
        edg->setPrev(lastEdge);
        lastEdge = edg;
    }

    return;
}

void TMesh::CreateEdges() {
    std::cout << "Creating edges... ";

    std::vector<TSide*> sides;
    std::vector<TSide*>::iterator it;

    TSide* sid = firstSide;

    FOR_SIDES(sid) {
        sid->setAux(false);
    }

    sid = firstSide;

    FOR_SIDES(sid) {
        if (sid->getAux()) {
            continue;
        }
        sid->setAux(true);
        sides.clear();
        sides.insert(sides.end(), sid);

        FOR_SIDE_NODES(sid->GetNNodes()) {

            FOR_NODE_SIDES(sid->getNode(iiSide), it) {
                if ((*it)->getAux()) {
                    continue;
                }
                if (sid->AreBBNeighbours((*it))) {
                    sides.insert(sides.end(), (*it));
                    (*it)->setAux(true);
                }
            }
        }
        TEdge* edg = new TEdge(sides);
        AddEdge(edg);
    }
    std::cout << n_edges << " edges created. OK.\n";
}

void TMesh::AddNeighbour(TNeighbour* ngh) {
    n_neighbours++;

    if (first_ngh == NULL) {
        first_ngh = ngh;
        last_ngh = ngh;
    } else {
        last_ngh->setNext(ngh);
        ngh->setPrev(last_ngh);
        last_ngh = ngh;
    }

    return;
}

void TMesh::CreateNeighboursBB() {
    std::cout << "  Creating BB neighbours... ";

    int count = 0;
    TEdge* edg;

    FOR_EDGES(edg, this) {
        if (edg->GetNSides() == 1) {
            continue;
        }

        TNeighbour* ngh = new TNeighbour(edg);
        AddNeighbour(ngh);
        count++;
    }

    std::cout << count << " neighbours created. OK.\n";
    return;
}

void TMesh::CreateNeighboursVB() {
    std::cout << "  Creating VB neighbours... ";

    int count = 0;
    TEdge* edg;

    FOR_EDGES(edg, this) {
        TElement* ele = NULL;
        std::vector<TElement *>::iterator ie;

        FOR_NODE_ELEMENTS(edg->getSide(0)->getNode(0), ie) {
            if (!edg->getSide(0)->AreVBNeighbours((*ie)))
                continue;
            else {
                ele = (*ie);
                break;
            }
        }
        if (ele != NULL) {
            TNeighbour* ngh = new TNeighbour(ele, edg->getSide(0));
            AddNeighbour(ngh);
            count++;
        }
    }
    std::cout << count << " neighbours created. OK.\n";
    return;
}

void TMesh::SideToNode() {
    std::cout << "Assigning sides to the nodes... ";

    TSide* sid = firstSide;

    FOR_SIDES(sid) {

        FOR_SIDE_NODES(sid->GetNNodes()) {
            TNode* nod = sid->getNode(iiSide);
            nod->addSide(sid);
        }
    }

    std::cout << "OK.\n";
    return;
}

void TMesh::CreateGeometry() {
    std::cout << "Creating geometry of the elements... ";

    TElement* ele;

    FOR_ELEMENTS(ele, this->getFirstElement()) {
        ele->CreateGeometry();
    }

    std::cout << "OK.\n";
}

void TMesh::TestElements() {
    std::cout << "Testing elements... ";

    TElement* ele;

    FOR_ELEMENTS(ele, this->getFirstElement()) {
        ele->Test();
    }

    std::cout << "OK.\n";
}
