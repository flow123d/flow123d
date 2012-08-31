#ifndef tetrahedronH
#define tetrahedronH

#include "point.h"
#include "triangle.h"

class TTetrahedron {
private:
    static int numberInstance;
    int id;

    TPoint* X1;
    TPoint* X2;
    TPoint* X3;
    TPoint* X4;

    TTriangle* T1;
    TTriangle* T2;
    TTriangle* T3;
    TTriangle* T4;

    TAbscissa* A1;
    TAbscissa* A2;
    TAbscissa* A3;
    TAbscissa* A4;
    TAbscissa* A5;
    TAbscissa* A6;

    double volume;

    int generateId();
    void ComputeVolume();

public:
    TTetrahedron(const TPoint&, const TPoint&, const TPoint&, const TPoint&);
    ~TTetrahedron();

    TTriangle GetTriangle(int) const;
    const TAbscissa &GetAbscissa(int) const;

    double GetMin(int) const;
    double GetMax(int) const;

    double GetVolume();

    static int getNumInstances() {
        return TTetrahedron::numberInstance;
    }
};
#endif
