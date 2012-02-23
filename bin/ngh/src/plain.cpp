#include "plain.h"
#include "mathfce.h"
#include "intersection.h"

using namespace mathfce;

int TPlain::numberInstance = 0;

int TPlain::generateId() {
    return TPlain::numberInstance++;
}

TPlain::TPlain() {
    id = generateId();

    U = new TVector(0, 0, 0);
    V = new TVector(0, 0, 0);
    N = new TVector(0, 0, 0);
    X = new TPoint(0, 0, 0);

    Compute();
}

TPlain::TPlain(const TVector &UU, const TVector &VV, const TPoint &XX) {
    id = generateId();

    U = new TVector(UU);
    V = new TVector(VV);
    N = new TVector(0, 0, 0);
    X = new TPoint(XX);

    Compute();
}

TPlain::TPlain(const TPoint &P1, const TPoint &P2, const TPoint &P3) {
    id = generateId();

    U = new TVector(P1, P2);
    V = new TVector(P1, P3);
    N = new TVector(0, 0, 0);
    X = new TPoint(P1);

    Compute();
}

TPlain::TPlain(const TPlain &p) {
    id = generateId();

    N = new TVector();
    U = new TVector();
    V = new TVector();
    X = new TPoint();

    a = p.a;
    b = p.b;
    c = p.c;
    d = p.d;

    *N = *p.N;
    *U = *p.U;
    *V = *p.V;
    *X = *p.X;
}

void TPlain::Compute() {
    *N = Cross(*U, *V);
    a = N->X1();
    b = N->X2();
    c = N->X3();
    d = -a * X->X() - b * X->Y() - c * X->Z();
    return;
}

TPlain::~TPlain() {
    delete U;
    delete V;
    delete N;
    delete X;
}

double TPlain::GetA() const {
    return a;
}

double TPlain::GetB() const {
    return b;
}

double TPlain::GetC() const {
    return c;
}

double TPlain::GetD() const {
    return d;
}

TVector TPlain::GetNormal() const {
    TVector tmp;
    tmp = *N;
    return tmp;
}

TVector TPlain::GetU() const {
    TVector tmp;
    tmp = *U;
    return tmp;
}

TVector TPlain::GetV() const {
    TVector tmp;
    tmp = *V;
    return tmp;
}

TPoint TPlain::GetPoint() const {
    TPoint tmp;
    tmp = *X;
    return tmp;
}

TPoint TPlain::GetPoint(double r, double s) const {
    TPoint tmp;
    tmp = r * *U + s * *V + *X;
    return tmp;
}

bool TPlain::Belong(const TPoint &P) const {
    if (IsZero(Distance(*this, P)))
        return true;
    return false;
}

TPlain & TPlain::operator =(const TPlain &p) {
    a = p.a;
    b = p.b;
    c = p.c;
    d = p.d;
    *(*this).U = *p.U;
    *(*this).V = *p.V;
    *(*this).N = *p.N;
    *(*this).X = *p.X;
    return *this;
}

void TPlain::SetPoints(const TPoint &P1, const TPoint &P2, const TPoint &P3) {
    *U = (TPoint) P2 - P1;
    *V = (TPoint) P3 - P1;
    *X = P1;
    Compute();
}

