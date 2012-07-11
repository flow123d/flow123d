
#include <iostream>

#include "Node.hh"

using namespace std;

Node::Node(const int label, const double x, const double y, const double z) {
    this->label = label;
    
    coor = new double[3];
    *(coor + 0) = x;
    *(coor + 1) = y;
    *(coor + 2) = z;
}

Node::~Node() {
//    cout << "Free node " << label << endl;
    delete coor;
}

int Node::getLabel() {
    return label;
}

void Node::setCoor(const int iCoor, const double value) {
    *(coor + iCoor) = value;
}

double Node::getCoor(const int iCoor) {
    return *(coor + iCoor);
}

void Node::setLevel(const double level) {
    this->level = level;
}

double Node::getLevel() {
    return level;
}

