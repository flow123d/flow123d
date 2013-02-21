#include "mesh/ngh/include/vertex.h"
#include "mesh/ngh/include/myvector.h"

int TVertex::numberInstance = 0;

int TVertex::generateId() {
    return TVertex::numberInstance++;
}

TVertex::TVertex(const TPoint& PP) {
    id = generateId();

    X = new TPoint(PP);
}

TVertex::~TVertex() {
    delete X;
}

TPoint TVertex::GetPoint() const {
    TPoint tmp;
    tmp = *X;
    return tmp;
}
