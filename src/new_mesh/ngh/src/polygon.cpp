#include <iostream>
#include <vector>

#include "new_mesh/ngh/include/polygon.h"
#include "new_mesh/ngh/include/vertex.h"
#include "new_mesh/ngh/include/mathfce.h"

using namespace mathfce;

int TPolygon::numberInstance = 0;

int TPolygon::generateId() {
    return TPolygon::numberInstance++;
}

TPolygon::TPolygon() {
    id = generateId();

    area_is_actual = false;
    center_is_actual = false;
    center.SetCoord(0, 0, 0);
    area = 0;
}

TPolygon::~TPolygon() {
    TPolygon *pol = this;
    std::vector<TVertex*>::iterator iv;

    FOR_POL_VERTECES(pol, iv) {
        delete (*iv);
    }
    verteces.clear();
}

std::ostream & operator <<(std::ostream& stream, const TPolygon& p) {
    std::vector<TVertex*>::iterator iv;
    int i = 0;

    // in the following, the explicit typecast is neccessary, since otherwise
    // the operators '=' and '!=' are not defined (I wander why constantness
    // of the Polygon instance can lead to this

    FOR_POL_VERTECES((TPolygon*) & p, iv) {
        stream << "V" << i << " = " << (*iv)->GetPoint() << "\n";
        i++;
    }
    return stream;
}

void TPolygon::Write() {
	std::vector<TVertex*>::iterator iv;
	printf("TPolygon::Write - size %d\n", verteces.size());
	FOR_POL_VERTECES(this, iv) {
		TPoint p = (*iv)->GetPoint();
		printf("%f %f %f\n", p.X(), p.Y(), p.Z());
	}
}

void TPolygon::Add(const TPoint& P) {
    std::vector<TVertex*>::iterator iv;

    FOR_POL_VERTECES(this, iv) {
        if ((*iv)->GetPoint() == P)
            return;
    }
    TVertex* V = new TVertex(P);
    int insertPos = InsertPosition(*V);
    area_is_actual = false;
    center_is_actual = false;
    verteces.insert(verteces.begin() + insertPos, V);
    return;
}

double TPolygon::GetArea() {
    if (!center_is_actual)
        ComputeCenter();
    printf("Center: %f %f %f\n", center.Get(1), center.Get(2), center.Get(3));
    getchar();
    if (area_is_actual)
        return area;
    ComputeArea();
    return area;
}

void TPolygon::ComputeCenter() {
    int i;
    double D[ 3 ] = {0, 0, 0};
    std::vector<TVertex*>::iterator iv;
    if (verteces.size() < 3) {
        center.SetCoord(D[ 0 ], D[ 1 ], D[ 2 ]);
        return;
    }

    FOR_POL_VERTECES(this, iv) {
        for (i = 0; i < 3; i++)
            D[ i ] += (*iv)->GetPoint().Get(i + 1);
    }
    for (i = 0; i < 3; i++)
        D[ i ] /= verteces.size();
    center.SetCoord(D[ 0 ], D[ 1 ], D[ 2 ]);
    center_is_actual = true;
    return;
}

TPoint TPolygon::GetCenter() {
    if (center_is_actual)
        return center;
    ComputeCenter();
    return center;
}

int TPolygon::InsertPosition(const TVertex& Vx) {
    TVertex *V1, *V2;
    TVector U1, U2;
    double d;
    int pos = 0;

    // if size < 2 new vertex is insert at the end
    if (verteces.size() < 2)
        return verteces.size();

    for (int i=0; i<verteces.size(); ++i) {
    	V1 = (*(verteces.begin() + i));
    	if (i == verteces.size() - 1)
    	    V2 = (*(verteces.begin()));
    	else
    	    V2 = (*(verteces.begin() + i + 1));
    	U1 = V2->GetPoint() - V1->GetPoint();
    	U2 = Vx.GetPoint() - V1->GetPoint();
    	d = U1.Get(1) * U2.Get(2) + U1.Get(2) * U2.Get(3) + U1.Get(3) * U2.Get(1)
    	        - U1.Get(3) * U2.Get(2) - U1.Get(2) * U2.Get(1) - U1.Get(1) * U2.Get(3);
    	if (verteces.size() == 2) {
    		if (d < 0.0) return 1;
    		else return 2;
    	}
    	if (d < 0) {
    		if (pos > 0) xprintf(Msg, " - TPolygon::InsertPosition - ERROR: only one vector product must be negative\n");
    		pos = i+1;
    	}
    }

    if (pos == 0) xprintf(Msg, " - TPolygon::InsertPosition - ERROR: no vector product is negative\n");
    return pos;
}

void TPolygon::ComputeArea() {
    std::vector<TVertex*>::iterator iv;

    TTriangle* T;
    TPoint P1, P2, P3;

    area = 0;
    area_is_actual = true;
    if (verteces.size() < 3)
        return;

    FOR_POL_VERTECES(this, iv) {
        P1 = center;
        P2 = (*iv)->GetPoint();
        if (iv == verteces.end() - 1)
            P3 = (*verteces.begin())->GetPoint();
        else
            P3 = (*(iv + 1))->GetPoint();
        T = new TTriangle(P1, P2, P3);
        area += T->GetArea();
    }
    return;
}

