
#include <math.h>

#include "Element.hh"
#include "ETetrahedron.hh"

using namespace std;

ETetrahedron::ETetrahedron(const int label, const short numTags, int* tags, Node** nodes) : Element(ELEMENT_TYPE_TETRAHEDRON, label, numTags, tags, nodes) {
    update();
}

ETetrahedron::~ETetrahedron() {
}

// =============================================================================

/**
 *                     |  1  x0  y0  z0  |
 * volume =  1/6 . det |  1  x1  y1  z1  |
 *                     |  1  x2  y2  z2  |
 *                     |  1  x3  y3  z3  |
 */

void ETetrahedron::update() {
    Node* node0 = getNode(0);
    double x0 = node0->getCoor(COOR_X);
    double y0 = node0->getCoor(COOR_Y);
    double z0 = node0->getCoor(COOR_Z);

    Node* node1 = getNode(1);
    double x1 = node1->getCoor(COOR_X);
    double y1 = node1->getCoor(COOR_Y);
    double z1 = node1->getCoor(COOR_Z);

    Node* node2 = getNode(2);
    double x2 = node2->getCoor(COOR_X);
    double y2 = node2->getCoor(COOR_Y);
    double z2 = node2->getCoor(COOR_Z);

    Node* node3 = getNode(3);
    double x3 = node3->getCoor(COOR_X);
    double y3 = node3->getCoor(COOR_Y);
    double z3 = node3->getCoor(COOR_Z);

    double value = 0.0;
    value += (x1 * y2 * z3) + (x2 * y3 * z1) + (x3 * y1 * z2) - ( x3 * y2 * z1 ) - (x2 * y1 * z3) - (x1 * y3 * z2);
    value -= (x0 * y2 * z3) + (x2 * y3 * z0) + (x3 * y0 * z2) - ( x3 * y2 * z0 ) - (x2 * y0 * z3) - (x0 * y3 * z2);
    value += (x0 * y1 * z3) + (x1 * y3 * z0) + (x3 * y0 * z1) - ( x3 * y1 * z0 ) - (x1 * y0 * z3) - (x0 * y3 * z1);
    value -= (x0 * y1 * z2) + (x1 * y2 * z0) + (x2 * y0 * z1) - ( x2 * y1 * z0 ) - (x1 * y0 * z2) - (x0 * y2 * z1);

    volume = fabs(value) / 6.0;
}

double ETetrahedron::getSize() {
    return volume;
}
