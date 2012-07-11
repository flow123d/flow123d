
#include <math.h>

#include "Element.hh"
#include "ETriangle.hh"

using namespace std;

ETriangle::ETriangle(const int label, const short numTags, int* tags, Node** nodes) : Element(ELEMENT_TYPE_TRIANGLE, label, numTags, tags, nodes) {
    update();
}

ETriangle::~ETriangle() {
}

// =============================================================================

/**
 *             |  1  x0  y0  |
 * sub0 =  det |  1  x1  y1  |
 *             |  1  x2  y2  |
 *
 *             |  1  y0  z0  |
 * sub1 =  det |  1  y1  z1  |
 *             |  1  y2  z2  |
 *
 *             |  1  z0  x0  |
 * sub2 =  det |  1  z1  x1  |
 *             |  1  z2  x2  |
 *
 * area =  1/2 . sqrt (  sub0^2 + sub1^2 + sub2^2  )
 *
 */

void ETriangle::update() {
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

    double sub0 = (x1 * y2) + (x2 * y0) + (x0 * y1) - (x1 * y0) - (x2 * y1) - (x0 * y2);
    double sub1 = (y1 * z2) + (y2 * z0) + (y0 * z1) - (y1 * z0) - (y2 * z1) - (y0 * z2);
    double sub2 = (z1 * x2) + (z2 * x0) + (z0 * x1) - (z1 * x0) - (z2 * x1) - (z0 * x2);

    area = 0.5 * sqrt((sub0 * sub0) + (sub1 * sub1) + (sub2 * sub2));
}

double ETriangle::getSize() {
    return area;
}
