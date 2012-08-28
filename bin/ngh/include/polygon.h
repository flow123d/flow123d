#ifndef polygonH
#define polygonH

#include "vertex.h"
#include "triangle.h"
#include <vector>

#define FOR_POL_VERTECES(i,j)     \
	for((j)=(i)->verteces.begin();(j)!=(i)->verteces.end();(j)++)

class TPolygon {
private:
    static int numberInstance;
    int id;

    double area;
    TPoint center;
    bool area_is_actual;
    bool center_is_actual;

    std::vector<TVertex*> verteces;

    int generateId();

    void ComputeCenter();
    void ComputeArea();

    /**
     * Find position of new vertex
     *
     * @param Vx New vertex
     * @return position of new vertex
     */
    int InsertPosition(const TVertex& Vx);

public:
    TPolygon();
    ~TPolygon();

    friend std::ostream & operator <<(std::ostream&, const TPolygon&);

    void Add(const TPoint&);
    void Write();
    double GetArea();
    TPoint GetCenter();

    static int getNumInstances() {
        return TPolygon::numberInstance;
    }
};

#endif

