
#include <memory>

#include "BoundaryData.hh"

using namespace std;

void BoundaryData::init(int id, int type, double value, int tag) {
    this->id = id;
    this->type = type;
    this->value = value;
    this->tag = tag;
}

BoundaryData::BoundaryData(int id, int type, double value, double sigma, int tag) {
    init(id, type, value, tag);
    this->sigma = sigma;
}

BoundaryData::BoundaryData(int id, int type, double value, int tag) {
    init(id, type, value, tag);
}

BoundaryData::~BoundaryData() {

}

int BoundaryData::getId() {
    return id;
}

int BoundaryData::getType() {
    return type;
}

double BoundaryData::getValue() {
    return value;
}

double BoundaryData::getSigma() {
    return sigma;
}

int BoundaryData::getTag() {
    return tag;
}
