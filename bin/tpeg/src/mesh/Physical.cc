
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "Physical.hh"

using namespace std;

Physical::Physical(const int dim, const int id, char* name) {
    this->dim = dim;
    this->id = id;

    this->name = (char *) malloc(sizeof (char) * (strlen(name) + 1));
    strcpy(this->name, name);
}

Physical::~Physical() {
    free(name);
}

int Physical::getDim() {
    return dim;
}

int Physical::getId() {
    return id;
}

char* Physical::getName() {
    char* name = (char *) malloc(sizeof (char) * (strlen(this->name) + 1));
    strcpy(name, this->name);
    return name;
}
