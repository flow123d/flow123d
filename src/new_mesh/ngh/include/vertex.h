#ifndef vertexH
#define vertexH

#include "point.h"

class TVertex {
private:
    static int numberInstance;
    int id;

    TPoint* X;

    int generateId();

public:
    TVertex(const TPoint&);
    ~TVertex();

    TPoint GetPoint() const;

    static int getNumInstances() {
        return TVertex::numberInstance;
    }
};

#endif

