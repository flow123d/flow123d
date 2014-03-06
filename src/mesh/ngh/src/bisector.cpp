#include <iostream>

#include "mesh/ngh/include/bisector.h"
#include "mesh/ngh/include/mathfce.h"
#include "mesh/ngh/include/intersection.h"
//#include "mesh/ngh/include/problem.h"

using namespace mathfce;

int TBisector::numberInstance = 0;

int TBisector::generateId() {
    return TBisector::numberInstance++;
}

TBisector::TBisector() {
    id = generateId();

    X0 = new TPoint(0, 0, 0);
    U = new TVector(0, 0, 0);
}

TBisector::TBisector(const TPoint &XX0, const TVector &UU) {
    id = generateId();

    X0 = new TPoint(XX0);
    U = new TVector(UU);
}

TBisector::TBisector(const TPoint &P0, const TPoint &P1) {
    id = generateId();

    X0 = new TPoint(P0);
    U = new TVector(P0, P1);
}

TBisector::TBisector(const  TBisector &x)
{
    id = generateId();

    X0 = new TPoint(*(x.X0));
    U = new TVector(*(x.U));
}

TBisector::~TBisector() {
    delete X0;
    delete U;
}

TBisector & TBisector::operator =(const TBisector &b) {
    //  U = new TVector();
    //  X0 = new TPoint();
    *(*this).U = *b.U;
    *(*this).X0 = *b.X0;

    return *this;
}

std::ostream & operator <<(std::ostream &stream, const TBisector &b) {
    stream << "U = (" << b.U->X1() << ", " << b.U->X2() << ", " << b.U->X3() << ")\n";
    stream << "X0 = [" << b.X0->X() << ", " << b.X0->Y() << ", " << b.X0->Z() << "]\n";

    return stream;
}

void TBisector::SetPoint(const TPoint &P) {
    *X0 = P;
}

void TBisector::SetVector(const TVector &UU) {
    *U = UU;
}

void TBisector::SetPoints(const TPoint &P0, const TPoint &P1) {
    *X0 = P0;
    *U = (TPoint) P1 - P0;
}

const TVector &TBisector::GetVector() const {

    return *(this->U);
}

const TPoint &TBisector::GetPoint() const
{
    return *(this->X0);
}

TPoint TBisector::GetPoint(double t) const {
    TPoint tmp;

    tmp.SetCoord(t * *U);
    tmp = tmp + *X0;

    return tmp;
}

void TBisector::GetParameter(const TPoint &P, double &t, bool &onBisector) const {
	t = (P.X() - X0->X()) / U->X1();
	onBisector = (fabs( (P.Y() - X0->Y()) / U->X2() - t ) < epsilon) & (fabs( (P.Z() - X0->Z()) / U->X3() - t ) < epsilon);
}

bool TBisector::Belong(const TPoint &P) const {
    if (IsZero(Distance(*this, P))) {
        return true;
    } else {
        return false;
    }
}
