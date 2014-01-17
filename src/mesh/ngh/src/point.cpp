#include <iostream>

#include "system/exc_common.hh"
#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/mathfce.h"

using namespace mathfce;

int TPoint::numberInstance = 0;

int TPoint::generateId() {
    return TPoint::numberInstance++;
}

TPoint::TPoint() {
    x = 0;
    y = 0;
    z = 0;

    id = generateId();
}

TPoint::TPoint(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;

    id = generateId();
}

TPoint::TPoint(const TPoint& point) {
    x = point.X();
    y = point.Y();
    z = point.Z();

    id = generateId();
}

TPoint::~TPoint() {
}

TPoint& TPoint::operator =(const TPoint& P) {
    x = P.x;
    y = P.y;
    z = P.z;

    return *this;
}
/*
bool TPoint::operator ==(TPoint* P) const {
    if (IsEqual(x, P->x) && IsEqual(y, P->y) && IsEqual(z, P->z)) {
        return true;
    }
    return false;
}
*/
bool TPoint::operator ==(const TPoint& P) const {
    if (IsEqual(x, P.x) && IsEqual(y, P.y) && IsEqual(z, P.z)) {
        return true;
    }
    return false;
}

TPoint& TPoint::operator =(const TVector& U) {
    x = U.X1();
    y = U.X2();
    z = U.X3();

    return *this;
}

TVector TPoint::operator -(const TPoint& P) const {
    return TVector( x - P.x, y - P.y, z - P.z);
}

/*
TPoint* TPoint::operator +(TPoint* P) const {
    TPoint* res = new TPoint();

    res->x = x + P->x;
    res->y = y + P->y;
    res->z = z + P->z;

    return res;
}
*/

TPoint TPoint::operator +(const TPoint& P) const {
    TPoint res;

    res.x = x + P.x;
    res.y = y + P.y;
    res.z = z + P.z;

    return res;
}

void TPoint::SetCoord(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

void TPoint::SetCoord(const TVector& U) {
    x = U.X1();
    y = U.X2();
    z = U.X3();
}

double TPoint::X() const {
    return x;
}

double TPoint::Y() const {
    return y;
}

double TPoint::Z() const {
    return z;
}

double TPoint::Get(int i) const {
    if (!(i >= 1 && i <= 3)) {
        THROW( ExcAssertMsg() << EI_Message( "Invalid specification of the element of the vector.") );
    }

    switch (i) {
        case 1: return x;
        case 2: return y;
        case 3: return z;
    }

    return 0.0;
}

std::ostream & operator <<(std::ostream& stream, const TPoint& P) {
    stream << "[" << P.x << " " << P.y << " " << P.z << "]";
    return stream;
}
