#include <cmath>

#include "config.h"
#include "tetrahedron.h"
#include "intersection.h"
#include "mathfce.h"
#include "system.h"

int TTetrahedron::numberInstance = 0;

int TTetrahedron::generateId() {
    return TTetrahedron::numberInstance++;
}

TTetrahedron::TTetrahedron(const TPoint& X1, const TPoint& X2,
        const TPoint& X3, const TPoint& X4) {
    id = generateId();

    this->X1 = new TPoint(X1);
    this->X2 = new TPoint(X2);
    this->X3 = new TPoint(X3);
    this->X4 = new TPoint(X4);

    T1 = new TTriangle(X2, X3, X4);
    T2 = new TTriangle(X1, X3, X4);
    T3 = new TTriangle(X1, X2, X4);
    T4 = new TTriangle(X1, X2, X3);

    A1 = new TAbscissa(X1, X2);
    A2 = new TAbscissa(X2, X3);
    A3 = new TAbscissa(X3, X1);
    A4 = new TAbscissa(X1, X4);
    A5 = new TAbscissa(X2, X4);
    A6 = new TAbscissa(X3, X4);

    ComputeVolume();
}

TTetrahedron::~TTetrahedron() {
    delete X1;
    delete X2;
    delete X3;
    delete X4;

    delete T1;
    delete T2;
    delete T3;
    delete T4;

    delete A1;
    delete A2;
    delete A3;
    delete A4;
    delete A5;
    delete A6;
}

TTriangle TTetrahedron::GetTriangle(int i) const {
    TTriangle tmp;
    switch (i) {
        case 1: tmp = *T1;
            break;
        case 2: tmp = *T2;
            break;
        case 3: tmp = *T3;
            break;
        case 4: tmp = *T4;
            break;
        default: mythrow((char*)"Unknown number of the triangle of the tetrahedron.", __LINE__, __FUNC__);
    }
    return tmp;
}

TAbscissa TTetrahedron::GetAbscissa(int i) const {
    TAbscissa tmp;
    switch (i) {
        case 1: tmp = *A1;
            break;
        case 2: tmp = *A2;
            break;
        case 3: tmp = *A3;
            break;
        case 4: tmp = *A4;
            break;
        case 5: tmp = *A5;
            break;
        case 6: tmp = *A6;
            break;
        default: mythrow((char*)"Unknown number of the triangle of the tetrahedron.", __LINE__, __FUNC__);
    }
    return tmp;
}

double TTetrahedron::GetMin(int i) const {
    double min = X1->Get(i);

    if (X2->Get(i) < min) {
        min = X2->Get(i);
    } else if (X3->Get(i) < min) {
        min = X3->Get(i);
    } else if (X4->Get(i) < min) {
        min = X4->Get(i);
    }

    return min;
}

double TTetrahedron::GetMax(int i) const {
    double max = X1->Get(i);

    if (X2->Get(i) > max) {
        max = X2->Get(i);
    } else if (X3->Get(i) > max) {
        max = X3->Get(i);
    } else if (X4->Get(i) > max) {
        max = X4->Get(i);
    }

    return max;
}

double TTetrahedron::GetVolume() {
    return volume;
}

void TTetrahedron::ComputeVolume() {
    double a[ 3 ][ 3 ];

    a[ 0 ][ 0 ] = X2->X() - X1->X();
    a[ 0 ][ 1 ] = X2->Y() - X1->Y();
    a[ 0 ][ 2 ] = X2->Z() - X1->Z();
    a[ 1 ][ 0 ] = X3->X() - X1->X();
    a[ 1 ][ 1 ] = X3->Y() - X1->Y();
    a[ 1 ][ 2 ] = X3->Z() - X1->Z();
    a[ 2 ][ 0 ] = X4->X() - X1->X();
    a[ 2 ][ 1 ] = X4->Y() - X1->Y();
    a[ 2 ][ 2 ] = X4->Z() - X1->Z();

    volume = fabs(Determinant3(a)) / 6.0;
}
