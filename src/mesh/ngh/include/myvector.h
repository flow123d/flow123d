#ifndef vectorH
#define vectorH

class TPoint;
#include "system.h"


class TVector {
protected:
    static int numberInstance;
    int id;

    double coors[ 3 ];
    double length;

    int generateId();

    void Compute();
    void CompLength();

public:
    TVector();
    TVector(double, double, double);
    TVector(TPoint, TPoint);
    TVector(const TVector &);
    ~TVector();

    double Length() const;
    void Get(double&, double&, double&) const;
    void Get(double*) const;
    double Get(int) const;
    void SetVector(double, double, double);
    bool IsZero();

    double X1() const;
    double X2() const;
    double X3() const;

    TVector & operator =(const TPoint&);
    TVector operator +(const TVector&);
    TVector operator +(const TPoint&);
    TVector operator -(const TVector&);
    friend TVector operator*(const TVector&, double);
    friend TVector operator*(double, const TVector&);
    bool operator ==(const TVector&);

    static int getNumInstances() {
        return TVector::numberInstance;
    }
};

TVector Cross(const TVector&, const TVector&);
double Dot(const TVector&, const TVector&);
bool AreParallel(const TVector&, const TVector&);
bool ArePerpendicular(const TVector&, const TVector&);

#endif
