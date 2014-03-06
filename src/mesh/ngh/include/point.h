#ifndef pointH
#define pointH

#include <iostream>
//#include "system.h"

#include "myvector.h"

class TVector;

class TPoint {
private:
    static int numberInstance;
    int id;

    double x;
    double y;
    double z;

    int generateId();

public:
    TPoint();
    TPoint(double, double, double);
    TPoint(const TPoint&);
    ~TPoint();

    TPoint * operator =(TPoint*);
    TPoint * operator +(TPoint*);
    bool operator ==(TPoint*);

    TPoint & operator =(const TPoint&);
    TPoint & operator =(const TVector&);
    TVector operator -(const TPoint&) const;
    TPoint operator +(const TPoint&) const;
    bool operator ==(const TPoint&) const;
    friend std::ostream & operator <<(std::ostream&, const TPoint&);

    void SetCoord(double, double, double);
    void SetCoord(const TVector&);

    double X() const;
    double Y() const;
    double Z() const;

    double Get(int) const;

    static int getNumInstances() {
        return TPoint::numberInstance;
    }
};

#endif
