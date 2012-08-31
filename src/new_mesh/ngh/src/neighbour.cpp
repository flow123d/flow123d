#include "new_mesh/ngh/include/intersectionLocal.h"
#include "new_mesh/ngh/include/neighbour.h"
#include "new_mesh/ngh/include/mesh.h"

int TNeighbour::numberInstance = 0;

int TNeighbour::generateId() {
    return TNeighbour::numberInstance++;
}


TNeighbour::TNeighbour(TEdge* edge)
: next(NULL), prev(NULL), intersection(NULL)
{
    id = generateId();
    type = nt_bb;

    n_sides = edge->GetNSides();
    side = new TSide*[n_sides];

    n_elements = edge->GetNSides();
    element = new TElement*[n_elements];

    this->edge = edge;

    int i;
    FOR_EDGE_SIDES(this->edge, i) {
        side[ i ] = this->edge->getSide(i);
        element[ i ] = side[ i ]->getElement();
    }
}

TNeighbour::TNeighbour(TElement* ele, TSide* sid)
: next(NULL), prev(NULL), intersection(NULL)
{
    id = generateId();
    type = nt_vb;

    n_sides = 2;
    side = new TSide*[n_sides];

    n_elements = 2;
    element = new TElement*[n_elements];

    edge = sid->getEdg();

    element[ 0 ] = ele;
    element[ 1 ] = sid->getElement();

    side[ 0 ] = NULL;
    side[ 1 ] = sid;
}

TNeighbour::TNeighbour(TElement* ele1, TElement* ele2, IntersectionLocal *intersec)
: prev(NULL), next(NULL), intersection(intersec)
{
    id = generateId();
    type = nt_vv;

    n_sides = 0;
    side = NULL;

    n_elements = 2;
    element = new TElement*[n_elements];

    edge = NULL;

    element[ 0 ] = ele1;
    element[ 1 ] = ele2;

}

TNeighbour::~TNeighbour() {
    if (side != NULL) {
        delete side;
        side = NULL;
    }
    if (element != NULL) {
        delete element;
        element = NULL;
    }
}

int TNeighbour::GetId() {
    return id;
}

int TNeighbour::GetNElements() {
    return n_elements;
}

TNghType TNeighbour::GetType() {
    return type;
}

/*
double TNeighbour::GetCoef() {
    return coef;
}*/

TSide* TNeighbour::getSide(int i) {
    return side[i];
}

TElement* TNeighbour::getElement(int i) {
    return element[i];
}

void TNeighbour::setPrev(TNeighbour* prev) {
    this->prev = prev;
}

TNeighbour* TNeighbour::getPrev() {
    return prev;
}

void TNeighbour::setNext(TNeighbour* next) {
    this->next = next;
}

TNeighbour* TNeighbour::getNext() {
    return next;
}
