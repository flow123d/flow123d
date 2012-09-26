#ifndef bisectorH
#define bisectorH

#include <iostream>

#include "point.h"
#include "myvector.h"

class TBisector {
protected:
    static int numberInstance;
    int id;

    TPoint* X0;
    TVector* U;

    int generateId();

public:
    TBisector();
    TBisector(const TPoint&, const TVector&);
    TBisector(const TPoint&, const TPoint&);
    TBisector(const  TBisector &);
    ~TBisector();

    TBisector & operator =(const TBisector&);
    friend std::ostream & operator <<(std::ostream&, const TBisector&);

    void SetPoints(const TPoint&, const TPoint&);
    bool Belong(const TPoint&) const;

    void SetPoint(const TPoint&);
    const TPoint &GetPoint() const;
    TPoint GetPoint(double) const;

    void SetVector(const TVector&);
    const TVector &GetVector() const;

    static int getNumInstances() {
        return TBisector::numberInstance;
    }
};

#endif
