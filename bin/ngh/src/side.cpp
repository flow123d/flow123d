#include <vector>

#include "side.h"

int TSide::numberInstance = 0;

int TSide::generateId() {
    return TSide::numberInstance++;
}

void TSide::Initialize() {
    shape = UNKNOWN;
    n_nodes = 0;

    prev = NULL;
    next = NULL;

    id = generateId();

    lnum = -1;
    ele = NULL;
    edg = NULL;
}

TSide::TSide(TElement* ele, int lnum, TNode* node1) {
    Initialize();

    shape = POINT;
    n_nodes = 1;

    node = new TNode*[ 1 ];
    node[ 0 ] = node1;

    this->lnum = lnum;
    this->ele = ele;
    ele->setSide(lnum, this);
}

TSide::TSide(TElement* ele, int lnum, TNode* node1, TNode* node2) {
    Initialize();

    shape = LINE;
    n_nodes = 2;

    node = new TNode*[ 2 ];
    node[ 0 ] = node1;
    node[ 1 ] = node2;

    this->lnum = lnum;
    this->ele = ele;
    ele->setSide(lnum, this);
}

TSide::TSide(TElement* ele, int lnum, TNode* node1, TNode* node2, TNode* node3) {
    Initialize();

    shape = TRIANGLE;
    n_nodes = 3;

    node = new TNode*[ 3 ];
    node[ 0 ] = node1;
    node[ 1 ] = node2;
    node[ 2 ] = node3;

    this->lnum = lnum;
    this->ele = ele;
    ele->setSide(lnum, this);
}

TSide::~TSide() {
    if (node != NULL) {
        delete node;
        node = NULL;
    }
}

int TSide::GetId() {
    return id;
}

bool TSide::AreBBNeighbours(TSide* sid) {
    if (shape != sid->shape) {
        return false;
    }
    for (int i = 0; i < n_nodes; i++) {
        bool test = false;
        for (int j = 0; j < sid->n_nodes; j++) {
            if (node[ i ] == sid->node[ j ]) {
                test = true;
            }
        }
        if (!test) {
            return false;
        }
    }

    TNode* nod = node[ 0 ];
    std::vector<TElement*>::iterator ei;

    FOR_NODE_ELEMENTS(nod, ei) {
        if (AreVBNeighbours((*ei))) {
            return false;
        }
    }
    return true;
}

bool TSide::AreVBNeighbours(TElement* ele) {
    if (shape != ele->GetType()) {
        return false;
    }
    for (int i = 0; i < n_nodes; i++) {
        bool test = false;
        for (int j = 0; j < ele->GetNNodes(); j++) {
            if (node[ i ] == ele->getNode(j)) {
                test = true;
            }
        }
        if (!test) {
            return false;
        }
    }
    return true;
}

int TSide::GetLnum() {
    return lnum;
}

int TSide::GetNNodes() {
    return n_nodes;
}

void TSide::setAux(bool aux) {
    this->aux = aux;
}

bool TSide::getAux() {
    return aux;
}

void TSide::setPrev(TSide* prev) {
    this->prev = prev;
}

void TSide::setNext(TSide* next) {
    this->next = next;
}

TSide* TSide::getNext() {
    return next;
}

void TSide::setEdg(TEdge* edg) {
    this->edg = edg;
}

TEdge* TSide::getEdg() {
    return edg;
}

TElement* TSide::getElement() {
    return ele;
}

TNode* TSide::getNode(int i) {
    return node[i];
}
