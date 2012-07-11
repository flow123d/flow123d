
#include <malloc.h>

#include"Element.hh"

using namespace std;

//                                  typ -  0  1  2   3  4   5   6   7   8   9  10  11  12  13  14 15
const short Element::NUMBER_OF_NODES[] = {-1, 2, 3, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1};

short Element::getNumNodes(const short elementType) {
    return NUMBER_OF_NODES[elementType];
}

Element::Element(const short elmType, const int label, const short numTags, int* tags, Node** nodes) {
    this->elementType = elmType;
    this->label = label;

    this->numTags = numTags;
    this->tags = (int*) malloc(numTags * sizeof (int));
    for (int iit = 0; iit < numTags; ++iit) {
        *(this->tags + iit) = *(tags + iit);
    }

    int numNodes = NUMBER_OF_NODES[elmType];
    this->nodes = (Node**) malloc(numNodes * sizeof (Node*));
    for (int iin = 0; iin < numNodes; ++iin) {
        *(this->nodes + iin) = *(nodes + iin);
    }

    swIsFictive = false;
}

Element::~Element() {
    free(nodes);
    free(tags);
}

int Element::getLabel() {
    return label;
}

short Element::getNumNodes() {
    return NUMBER_OF_NODES[elementType];
}

Node* Element::getNode(const short iin) {
    return *(nodes + iin);
}

short Element::getElementType() {
    return elementType;
}

short Element::getNumTags() {
    return numTags;
}

int Element::getTag(const short iTag) {
    return tags[iTag];
}

void Element::setFictive(const bool swIsFictive) {
    this->swIsFictive = swIsFictive;
}

bool Element::isFictive() {
    return swIsFictive;
}
