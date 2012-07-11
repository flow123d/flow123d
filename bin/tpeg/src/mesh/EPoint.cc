
#include "Element.hh"
#include "EPoint.hh"

using namespace std;

EPoint::EPoint(const int label, const short numTags, int* tags, Node** nodes) : Element(ELEMENT_TYPE_POINT, label, numTags, tags, nodes) {
    size = 0.0;
}

EPoint::~EPoint() {
}

// =============================================================================

void EPoint::update() {
}

double EPoint::getSize() {
    return size;
}
