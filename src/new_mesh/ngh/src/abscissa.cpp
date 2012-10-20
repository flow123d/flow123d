#include <cmath>

#include "new_mesh/ngh/include/abscissa.h"
#include "stdio.h"

int TAbscissa::numberInstance = 0;

int TAbscissa::generateId() {
    return TAbscissa::numberInstance++;
}

TAbscissa::TAbscissa() {
    id = generateId();
}

TAbscissa::TAbscissa(const TPoint& PP0, const TPoint& PP1) : TBisector(PP0, PP1) {
    id = generateId();

//    P0 = new TPoint(PP0);
//    P1 = new TPoint(PP1);

    //    *(*this)P0 = *PP0;
    //    *P1 = PP1;

    ComputeLength();
}

TAbscissa & TAbscissa::operator =(const TAbscissa& a) {
    //  U = new TVector();
    //  X0 = new TPoint();
    //  P0 = new TPoint();
    //  P1 = new TPoint();

    *U = *(a.U);
    *X0 = *(a.X0);

//    *(*this).P0 = *a.P0;
//    *(*this).P1 = *a.P1;

    length = a.length;

    return *this;
}

TAbscissa::~TAbscissa() {
//    delete P0;
//    delete P1;
}

void TAbscissa::SetPoints(const TPoint& PP0, const TPoint& PP1) {
//    *P0 = PP0;
//    *P1 = PP1;

    TBisector::SetPoints(PP0, PP1);
    ComputeLength();
}

void TAbscissa::ComputeLength() {

    length = U->Length();
}

double TAbscissa::Length() {
    return length;
}

BoundingBox &TAbscissa::get_bounding_box() {
	arma::vec3 minCoor;
	arma::vec3 maxCoor;

	for (int i=0; i<3; ++i) {
		minCoor(i) = GetMin(i+1);
		maxCoor(i) = GetMax(i+1);
	}

	boundingBox.set_bounds(minCoor, maxCoor);

	return boundingBox;
}

double TAbscissa::GetMin(int x) const {
    if (U == NULL ) printf("U is null!");
    if (U->Get(x) < 0 ) {
        return U->Get(x)+X0->Get(x);
    } else {
        return X0->Get(x);
    }
}

double TAbscissa::GetMax(int x) const {
    if (U->Get(x) > 0 ) {
        return U->Get(x)+X0->Get(x);
    } else {
        return X0->Get(x);
    }
}
