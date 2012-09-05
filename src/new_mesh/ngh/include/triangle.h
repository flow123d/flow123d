#ifndef triangleH
#define triangleH

#include "point.h"
#include "plain.h"
#include "abscissa.h"
#include "new_mesh/bounding_box.hh"

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

    BoundingBox* boundingBox;

    double area;

    int generateId();

    void ComputeArea();
    void compute_bounding_box();

public:
    TTriangle();
    TTriangle(const TTriangle&);
    TTriangle(const TPoint&, const TPoint&, const TPoint&);
    ~TTriangle();

    TTriangle & operator =(const TTriangle &t);

    const TPlain &GetPlain() const;
    const TAbscissa &GetAbscissa(int) const;
    const TPoint &GetPoint(int) const;

    void SetPoints(const TPoint&, const TPoint&, const TPoint&);

    double GetMin(int) const;
    double GetMax(int) const;

    double GetArea();
    BoundingBox* get_bounding_box();

    bool IsInner(const TPoint&) const;

    static int getNumInstances() {
        return TTriangle::numberInstance;
    }
};

#endif

