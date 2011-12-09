#include <iostream>
#include <vector>

#include "polygon.h"
#include "vertex.h"
#include "mathfce.h"

using namespace mathfce;

int TPolygon::numberInstance = 0;

int TPolygon::generateId() {
    return TPolygon::numberInstance++;
}

TPolygon::TPolygon() {
    id = generateId();

    area_is_actual = false;
    center_is_actual = false;
    vertex_order_is_actual = false;
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

void TPolygon::Add(const TPoint& P) {
    std::vector<TVertex*>::iterator iv;

    FOR_POL_VERTECES(this, iv) {
        if ((*iv)->GetPoint() == P)
            return;
    }
    TVertex* V = new TVertex(P);
    area_is_actual = false;
    center_is_actual = false;
    vertex_order_is_actual = false;
    verteces.insert(verteces.end(), V);
    return;
}

double TPolygon::GetArea() {
    if (!center_is_actual)
        ComputeCenter();
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

void TPolygon::TestVertexOrder() {
    std::vector<TVertex*>::iterator iv;
    TVertex *V1, *V2, *V3;
    TVector U1, U2;
    double d;
    if (vertex_order_is_actual)
        return;
    if (verteces.size() < 3)
        return;

    FOR_POL_VERTECES(this, iv) {
        if (iv == verteces.end() - 1)
            break;
        V1 = (*iv);
        V2 = (*(iv + 1));
        if (iv + 1 == verteces.end() - 1)
            V3 = (*(verteces.begin()));
        else
            V3 = (*(iv + 2));
        U1 = V2->GetPoint() - V1->GetPoint();
        U2 = V3->GetPoint() - V1->GetPoint();
        d = U1.Get(1) * U2.Get(2) + U1.Get(2) * U2.Get(3) + U1.Get(3) * U2.Get(1)
                - U1.Get(3) * U2.Get(2) - U1.Get(2) * U2.Get(1) - U1.Get(1) * U2.Get(3);
        if (d < 0) {
            if (V3 != (*verteces.begin()))
                iter_swap(iv + 1, iv + 2);
            else
                iter_swap(iv, iv + 1);
            iv = verteces.begin() - 1;
        }
    }
    vertex_order_is_actual = true;
}

void TPolygon::ComputeArea() {
    std::vector<TVertex*>::iterator iv;

    TTriangle T;

    TPoint P1, P2, P3;
    area = 0;
    if (verteces.size() < 3)
        return;
    TestVertexOrder();

    FOR_POL_VERTECES(this, iv) {
        P1 = center;
        P2 = (*iv)->GetPoint();
        if (iv == verteces.end() - 1)
            P3 = (*verteces.begin())->GetPoint();
        else
            P3 = (*(iv + 1))->GetPoint();
        T.SetPoints(P1, P2, P3);
        area += T.GetArea();
    }
    return;
}

