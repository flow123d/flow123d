
#include "Boundary.hh"

using namespace std;

void Boundary::init(int id, int type, double value, int where, int elmId, int sidId, int tag) {
    this->id = id;
    this->type = type;
    this->value = value;
    this->where = where;
    this->elmId = elmId;
    this->sidId = sidId;
    this->tag = tag;
}

Boundary::Boundary(int id, int type, double value, double sigma, int where, int elmId, int sidId, int tag) {
    init(id, type, value, where, elmId, sidId, tag);
    this->sigma = sigma;
}

Boundary::Boundary(int id, int type, double value, int where, int elmId, int sidId, int tag) {
    init(id, type, value, where, elmId, sidId, tag);
}

Boundary::~Boundary() {

}

int Boundary::getId() {
    return id;
}

int Boundary::getType() {
    return type;
}

double Boundary::getValue() {
    return value;
}

double Boundary::getSigma() {
    return sigma;
}

int Boundary::getWhere() {
    return where;
}

int Boundary::getElmId() {
    return elmId;
}

int Boundary::getSidId() {
    return sidId;
}

int Boundary::getTag() {
    return tag;
}
