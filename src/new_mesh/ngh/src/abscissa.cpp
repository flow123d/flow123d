#include <cmath>

#include "new_mesh/ngh/include/abscissa.h"

int TAbscissa::numberInstance = 0;

int TAbscissa::generateId() {
    return TAbscissa::numberInstance++;
}

TAbscissa::TAbscissa() {
    id = generateId();
}

TAbscissa::TAbscissa(const TPoint& PP0, const TPoint& PP1) : TBisector(PP0, PP1) {
    id = generateId();

    P0 = new TPoint(PP0);
    P1 = new TPoint(PP1);

    //    *(*this)P0 = *PP0;
    //    *P1 = PP1;

    ComputeLength();
}

TAbscissa & TAbscissa::operator =(const TAbscissa& a) {
    //  U = new TVector();
    //  X0 = new TPoint();
    //  P0 = new TPoint();
    //  P1 = new TPoint();

    *(*this).U = *a.U;

    *(*this).X0 = *a.X0;

    *(*this).P0 = *a.P0;
    *(*this).P1 = *a.P1;

    length = a.length;

    return *this;
}

TAbscissa::~TAbscissa() {
    delete P0;
    delete P1;
}

void TAbscissa::SetPoints(const TPoint& PP0, const TPoint& PP1) {
    *P0 = PP0;
    *P1 = PP1;

    TBisector::SetPoints(PP0, PP1);
}

void TAbscissa::ComputeLength() {
    double dx = (P1->X() - P0->X())*(P1->X() - P0->X());
    double dy = (P1->Y() - P0->Y())*(P1->Y() - P0->Y());
    double dz = (P1->Z() - P0->Z())*(P1->Z() - P0->Z());

    length = sqrt(dx * dx + dy * dy + dz * dz);
}

double TAbscissa::Length() {
    return length;
}

double TAbscissa::GetMin(int x) const {
    if (P0->Get(x) < P1->Get(x)) {
        return P0->Get(x);
    } else {
        return P1->Get(x);
    }
}

double TAbscissa::GetMax(int x) const {
    if (P0->Get(x) > P1->Get(x)) {
        return P0->Get(x);
    } else {
        return P1->Get(x);
    }
}
