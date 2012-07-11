
#include <math.h>
#include <iostream>

#include "Element.hh"
#include "ELine.hh"

using namespace std;

ELine::ELine(const int label, const short numTags, int* tags, Node** nodes) : Element(ELEMENT_TYPE_LINE, label, numTags, tags, nodes) {
    update();
}

ELine::~ELine() {
}

// =============================================================================

void ELine::update() {
    Node* node0 = getNode(0);
    double x0 = node0->getCoor(COOR_X);
    double y0 = node0->getCoor(COOR_Y);
    double z0 = node0->getCoor(COOR_Z);

    Node* node1 = getNode(1);
    double x1 = node1->getCoor(COOR_X);
    double y1 = node1->getCoor(COOR_Y);
    double z1 = node1->getCoor(COOR_Z);

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;

    length = sqrt(dx * dx + dy * dy + dz * dz);
}

double ELine::getSize() {
    return length;
}
