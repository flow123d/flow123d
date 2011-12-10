#include "neighbour.h"
#include "mesh.h"

int TNeighbour::numberInstance = 0;

int TNeighbour::generateId() {
    return TNeighbour::numberInstance++;
}

void TNeighbour::Initialize() {
    id = generateId();

    type = (TNghType) 0;

    n_sides = 0;
    side = NULL;

    n_elements = 0;
    element = NULL;

    edge = NULL;

    prev = NULL;
    next = NULL;

    coef = 0.0;
}

TNeighbour::TNeighbour(TEdge* edge) {
    Initialize();

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

    return;
}

TNeighbour::TNeighbour(TElement* ele, TSide* sid) {
    Initialize();

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

    coef = 1.0;

    return;
}

TNeighbour::TNeighbour(TElement* ele1, TElement* ele2, double coef) {
    Initialize();

    type = nt_vv;

    n_sides = 0;
    side = NULL;

    n_elements = 2;
    element = new TElement*[n_elements];

    edge = NULL;

    element[ 0 ] = ele1;
    element[ 1 ] = ele2;

    this->coef = coef;

    return;
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

double TNeighbour::GetCoef() {
    return coef;
}

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
