#include <cmath>
#include "new_mesh/ngh/include/point.h"
#include "new_mesh/ngh/include/myvector.h"
//#include <math.h>
#include "new_mesh/ngh/include/mathfce.h"

using namespace mathfce;

double epsilon;

int TVector::numberInstance = 0;

int TVector::generateId() {
    return TVector::numberInstance++;
}

TVector::TVector() {
    id = generateId();

    coors[ 0 ] = 0;
    coors[ 1 ] = 0;
    coors[ 2 ] = 0;

    Compute();
}

TVector::TVector(double xx1, double xx2, double xx3) {
    id = generateId();

    coors[ 0 ] = xx1;
    coors[ 1 ] = xx2;
    coors[ 2 ] = xx3;

    Compute();
}

TVector::TVector(TPoint P1, TPoint P2) {
    id = generateId();

    coors[ 0 ] = P2.X() - P1.X();
    coors[ 1 ] = P2.Y() - P1.Y();
    coors[ 2 ] = P2.Z() - P1.Z();

    Compute();
}

TVector::TVector(const TVector &x)
{
    id = generateId();

    coors[0]=x.coors[0];
    coors[1]=x.coors[1];
    coors[2]=x.coors[2];

    Compute();
}

TVector::~TVector() {
    ;
}

void TVector::Compute() {
    CompLength();
}

void TVector::CompLength() {
    length = sqrt(pow(this->coors[ 0 ], 2.0) + pow(this->coors[ 1 ], 2.0) + pow(this->coors[ 2 ], 2.0));
}

TVector& TVector::operator =(const TPoint &P) {
    this->coors[ 0 ] = P.X();
    this->coors[ 1 ] = P.Y();
    this->coors[ 2 ] = P.Z();

    Compute();

    return *this;
}

TVector TVector::operator +(const TVector &V) {
    TVector res(0, 0, 0);

    for (int i = 0; i < 3; i++) {
        res.coors[ i ] = coors[i] + V.coors[ i ];
    }

    res.Compute();

    return res;
}

TVector TVector::operator +(const TPoint &P) {
    TVector res(0, 0, 0);

    for (int i = 0; i < 3; i++) {
        res.coors[ i ] = coors[i] + P.Get(i + 1);
    }

    res.Compute();

    return res;
}

TVector TVector::operator -(const TVector &V) {
    TVector res(0, 0, 0);

    for (int i = 0; i < 3; i++) {
        res.coors[ i ] = coors[i] - V.coors[ i ];
    }

    res.Compute();

    return res;
}

void TVector::SetVector(double xx1, double xx2, double xx3) {
    coors[ 0 ] = xx1;
    coors[ 1 ] = xx2;
    coors[ 2 ] = xx3;
    Compute();
}

TVector operator *(const TVector &U, double x) {
    TVector tmp;
    tmp.SetVector(U.X1() * x, U.X2() * x, U.X3() * x);
    return tmp;
}

TVector operator *(double x, const TVector &U) {
    TVector tmp;
    tmp.SetVector(U.X1() * x, U.X2() * x, U.X3() * x);
    return tmp;
}

bool TVector::IsZero() {
    if (!mathfce::IsZero(length)) {
        return false;
    } else {
        return true;
    }
}

double TVector::Length() const {
    return length;
}

void TVector::Get(double &xx1, double &xx2, double &xx3) const {
    xx1 = coors[ 0 ];
    xx2 = coors[ 1 ];
    xx3 = coors[ 2 ];
    return;
}

void TVector::Get(double *U) const {
    for (int i = 0; i < 3; i++) {
        U[ i ] = coors[ i ];
    }

    return;
}

double TVector::Get(int i) const {
    return coors[ i - 1 ];
}

TVector Cross(const TVector &U, const TVector &V) {
    double x1, x2, x3;
    double u1, u2, u3;
    double v1, v2, v3;

    U.Get(u1, u2, u3);
    V.Get(v1, v2, v3);

    x1 = u2 * v3 - u3 * v2;
    x2 = u3 * v1 - u1 * v3;
    x3 = u1 * v2 - u2 * v1;

    TVector X(x1, x2, x3);

    return X;
}

double Dot(const TVector &U, const TVector &V) {
    double u1, u2, u3;
    double v1, v2, v3;

    U.Get(u1, u2, u3);
    V.Get(v1, v2, v3);

    double product = u1 * v1 + u2 * v2 + u3 * v3;

    return product;
}

bool TVector::operator ==(const TVector &U) {
    for (int i = 0; i < 3; i++) {
        if ((fabs(coors[ i ]) - fabs(U.coors[ i ])) > epsilon) {
            return false;
        }
    }
    return true;
}

double TVector::X1() const {
    return coors[ 0 ];
}

double TVector::X2() const {
    return coors[ 1 ];
}

double TVector::X3() const {
    return coors[ 2 ];
}

bool AreParallel(const TVector &U, const TVector &V) {
    // if vector W is zero, then these two bisectors are parallel or same
    /*  TVector W;
      W = (TVector)U - Dot(U, V) / Dot(V, V) * V;
      if (W.IsZero())
        return true;
      else
        return false;*/
    TVector W;
    W = Cross(U, V);
    if ((W.Length() < epsilon * 1 * U.Length())
            || (W.Length() < epsilon * 1 * V.Length())) {
        return true;
    } else {
        return false;
    }
}

bool ArePerpendicular(const TVector &U, const TVector &V) {
    double x = Dot(U, V) / (U.Length() * V.Length());
    return IsZero(x);
}
