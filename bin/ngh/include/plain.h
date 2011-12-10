#ifndef plainH
#define plainH

#include "myvector.h"
#include "point.h"

class TPlain {
private:
    static int numberInstance;
    int id;

    TVector* U;
    TVector* V;
    TVector* N;

    double a;
    double b;
    double c;
    double d;

    TPoint* X;

    int generateId();

    void Compute();

public:
    TPlain();
    TPlain(const TVector&, const TVector&, const TPoint&);
    TPlain(const TPoint&, const TPoint&, const TPoint&);
    TPlain(const TPlain&);
    ~TPlain();

    TPoint GetPoint() const;
    TPoint GetPoint(double, double) const;
    TVector GetNormal() const;
    TVector GetU() const;
    TVector GetV() const;
    double GetA() const;
    double GetB() const;
    double GetC() const;
    double GetD() const;
    bool Belong(const TPoint&) const;
    void SetPoints(const TPoint&, const TPoint&, const TPoint&);

    TPlain & operator =(const TPlain&);

    static int getNumInstances() {
        return TPlain::numberInstance;
    }
};

#endif
