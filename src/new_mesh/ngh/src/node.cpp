#include "new_mesh/ngh/include/config.h"
#include "new_mesh/ngh/include/node.h"
#include "new_mesh/ngh/include/system.h"

int TNode::numberInstance = 0;

int TNode::generateId() {
    return TNode::numberInstance++;
}

void TNode::Initialize() {
    id = generateId();

    label = -1;
    next = NULL;
    prev = NULL;
}

void TNode::setLabel(int label) {
    if (label < 0) {
        mythrow((char*) "The label of the node has to be >= 0.", __LINE__, __FUNC__);
    }
    this->label = label;
}

void TNode::ParseLine(char *line) {
    setLabel(atoi(strtok(line, " \t")));

    double x = atof(strtok(NULL, " \t"));
    double y = atof(strtok(NULL, " \t"));
    double z = atof(strtok(NULL, " \t"));

    SetCoord(x, y, z);
}

TNode::TNode() {
    Initialize();
}

TNode::TNode(double x, double y, double z) {
    Initialize();
    SetCoord(x, y, z);
}

TNode::TNode(int label, double x, double y, double z) {
    Initialize();
    this->label = label;
    SetCoord(x, y, z);
}

TNode::~TNode() {
}

int TNode::getLabel() {
    return label;
}

void TNode::setPrev(TNode* prev) {
    this->prev = prev;
}

void TNode::setNext(TNode* next) {
    this->next = next;
}

TNode* TNode::getNext() {
    return next;
}

std::vector<TSide*>::iterator TNode::getSideBegin() {
    return sides.begin();
}

std::vector<TSide*>::iterator TNode::getSideEnd() {
    return sides.end();
}

void TNode::addSide(TSide* side) {
    sides.insert(sides.end(), side);
}

std::vector<TElement*>::iterator TNode::getElementBegin() {
    return elements.begin();
}

std::vector<TElement*>::iterator TNode::getElementEnd() {
    return elements.end();
}

void TNode::addElement(TElement* ele) {
    elements.insert(elements.end(), ele);
}
