#include "new_mesh/ngh/include/triangle.h"
#include "new_mesh/ngh/include/config.h"
#include "new_mesh/ngh/include/system.h"


int TTriangle::numberInstance = 0;

int TTriangle::generateId() {
    return TTriangle::numberInstance++;
}

TTriangle::TTriangle() {
    id = generateId();

    area = 0;
}

TTriangle::TTriangle(const TTriangle& T) {
    id = generateId();

    X1 = new TPoint(*T.X1);
    X2 = new TPoint(*T.X2);
    X3 = new TPoint(*T.X3);

    A1 = new TAbscissa(*T.A1);
    A2 = new TAbscissa(*T.A2);
    A3 = new TAbscissa(*T.A3);

    pl = new TPlain(*T.pl);

    area = T.area;
    compute_bounding_box();
}

TTriangle::TTriangle(const TPoint& P1, const TPoint& P2, const TPoint& P3) {
    id = generateId();

    X1 = new TPoint(P1);
    X2 = new TPoint(P2);
    X3 = new TPoint(P3);

    A1 = new TAbscissa(P1, P2);
    A2 = new TAbscissa(P2, P3);
    A3 = new TAbscissa(P3, P1);

    pl = new TPlain(P1, P2, P3);

    ComputeArea();
    compute_bounding_box();
}

TTriangle::~TTriangle() {
    if (X1 != NULL) {
        delete X1;
    }
    if (X2 != NULL) {
        delete X2;
    }
    if (X3 != NULL) {
        delete X3;
    }

    if (A1 != NULL) {
        delete A1;
    }
    if (A2 != NULL) {
        delete A2;
    }
    if (A3 != NULL) {
        delete A3;
    }

    if (pl != NULL) {
        delete pl;
    }
}

const TPlain &TTriangle::GetPlain() const {
    return *pl;
}

const TAbscissa &TTriangle::GetAbscissa(const int i) const {
    switch (i) {
        case 1: return *A1;
            break;
        case 2: return *A2;
            break;
        case 3: return *A3;
            break;
        default:
            mythrow((char*) "Unknown number of the abscissa of the triangle.",  __LINE__, __FUNC__);

    }
}

const TPoint &TTriangle::GetPoint(int i) const {
    switch (i) {
        case 1: return *X1;
            break;
        case 2: return *X2;
            break;
        case 3: return *X3;
            break;
        default: mythrow((char*) "Unknown number of the point of the triangle.", __LINE__, __FUNC__);
    }
}

void TTriangle::SetPoints(const TPoint& P1, const TPoint& P2, const TPoint& P3) {
    *X1 = P1;
    *X2 = P2;
    *X3 = P3;

    A1->SetPoints(P1, P2);
    A2->SetPoints(P2, P3);
    A3->SetPoints(P3, P1);

    pl->SetPoints(P1, P2, P3);

    ComputeArea();
    compute_bounding_box();
}

void TTriangle::ComputeArea() {
    TVector N = Cross(A1->GetVector(), A2->GetVector());
    area = 0.5 * N.Length();
}

void TTriangle::compute_bounding_box() {
	arma::vec3 minCoor;
	arma::vec3 maxCoor;

	for (int i=0; i<3; ++i) {
		minCoor(i) = GetMin(i+1);
		maxCoor(i) = GetMax(i+1);
	}

	boundingBox = new BoundingBox(minCoor, maxCoor);
}

double TTriangle::GetArea() {
    return area;
}

const BoundingBox &TTriangle::get_bounding_box() const {
	return *boundingBox;
}

double TTriangle::GetMin(int i) const {
    double min = X1->Get(i);

    if (X2->Get(i) < min) {
        min = X2->Get(i);
    }
    if (X3->Get(i) < min) {
        min = X3->Get(i);
    }

    return min;
}

double TTriangle::GetMax(int i) const {
    double max = X1->Get(i);

    if (X2->Get(i) > max) {
        max = X2->Get(i);
    }
    if (X3->Get(i) > max) {
        max = X3->Get(i);
    }

    return max;
}

TTriangle & TTriangle::operator =(const TTriangle& t) {
    area = t.area;

    *(*this).A1 = *t.A1;
    *(*this).A2 = *t.A2;
    *(*this).A3 = *t.A3;
    *(*this).pl = *t.pl;
    *(*this).X1 = *t.X1;
    *(*this).X2 = *t.X2;
    *(*this).X3 = *t.X3;

    return *this;
}

bool TTriangle::IsInner(const TPoint& P) const {
    TVector N1, N2, U1(*X1, *X2), U2(*X2, *X3), U3(*X3, *X1);
    TVector Up1(*X1, P), Up2(*X2, P), Up3(*X3, P);

    N1 = Cross(Up1, U1);
    N2 = Cross(U1, U3);
    if (Dot(N1, N2) < 0) {
        return false;
    }

    N1 = Cross(Up2, U2);
    N2 = Cross(U2, U1);
    if (Dot(N1, N2) < 0) {
        return false;
    }

    N1 = Cross(Up3, U3);
    N2 = Cross(U3, U2);
    if (Dot(N1, N2) < 0) {
        return false;
    }

    return true;
}
