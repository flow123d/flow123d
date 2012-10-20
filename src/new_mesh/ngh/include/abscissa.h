#ifndef abscissaH
#define abscissaH

#include "bisector.h"
#include "point.h"
#include "new_mesh/bounding_box.hh"

class TAbscissa : public TBisector {
private:
    static int numberInstance;
    int id;

//    TPoint* P0;
//    TPoint* P1;

    BoundingBox boundingBox;

    double length;

    int generateId();
    void ComputeLength();

public:
    TAbscissa();
//    TAbscissa(double, double);
    TAbscissa(const TPoint&, const TPoint&);
    ~TAbscissa();

    TAbscissa & operator =(const TAbscissa&);

    double Length();
    BoundingBox &get_bounding_box();

    void SetPoints(const TPoint&, const TPoint&);

    double GetMin(int) const;
    double GetMax(int) const;

    static int getNumInstances() {
        return TAbscissa::numberInstance;
    }
};

#endif
