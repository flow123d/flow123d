#ifndef triangleH
#define triangleH

#include "point.h"
#include "plain.h"
#include "abscissa.h"

class TTriangle {
private:
    static int numberInstance;
    int id;

    TPoint* X1;
    TPoint* X2;
    TPoint* X3;

    TAbscissa* A1;
    TAbscissa* A2;
    TAbscissa* A3;

    TPlain* pl;

    double area;

    int generateId();

    void ComputeArea();

public:
    TTriangle();
    TTriangle(const TTriangle&);
    TTriangle(const TPoint&, const TPoint&, const TPoint&);
    ~TTriangle();

    TTriangle & operator =(const TTriangle &t);

    TPlain GetPlain() const;
    TAbscissa GetAbscissa(int) const;
    TPoint GetPoint(int) const;

    void SetPoints(const TPoint&, const TPoint&, const TPoint&);

    double GetMin(int) const;
    double GetMax(int) const;

    double GetArea();

    bool IsInner(const TPoint&) const;

    static int getNumInstances() {
        return TTriangle::numberInstance;
    }
};

#endif

