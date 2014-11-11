#ifndef abscissaH
#define abscissaH

#include "bisector.h"
#include "point.h"
#include "mesh/bounding_box.hh"

class TAbscissa : public TBisector {
private:
    static int numberInstance;
    int id;

    BoundingBox boundingBox;

    double length;

    int generateId();
    void ComputeLength();

public:
    TAbscissa();
    TAbscissa(const TPoint&, const TPoint&);
    TAbscissa(const Element&);
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
