#include <vector>

#include "edge.h"
#include "side.h"

int TEdge::numberInstance = 0;

TEdge::TEdge(std::vector<TSide *> sides) {
    id = TEdge::numberInstance++;

    prev = NULL;
    next = NULL;

    n_sides = sides.size();
    sid = new int[sides.size()];
    side = new TSide*[sides.size()];

    std::vector<TSide *>::iterator it;

    int i = 0;
    for (it = sides.begin(); it != sides.end(); it++, i++) {
        sid[ i ] = (*it)->GetLnum();
        side[ i ] = (*it);
        (*it)->setEdg(this);
    }
}

TEdge::~TEdge() {
    if (sid != NULL) {
        delete sid;
        sid = NULL;
    }
    if (side != NULL) {
        delete side;
        side = NULL;
    }
}

int TEdge::GetId() {
    return id;
}

int TEdge::GetNSides() {
    return n_sides;
}

int TEdge::getSId(int i) {
    return sid[i];
}

void TEdge::setPrev(TEdge* prev) {
    this->prev = prev;
}

TEdge* TEdge::getPrev() {
    return prev;
}

void TEdge::setNext(TEdge* next) {
    this->next = next;
}

TEdge* TEdge::getNext() {
    return next;
}

TSide* TEdge::getSide(int i) {
    return side[i];
}
