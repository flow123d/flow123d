#include "new_mesh/ngh/include/intersection.h"
#include "new_mesh/ngh/include/matrix.h"
#include "new_mesh/ngh/include/mathfce.h"
#include <cmath>
#include <iostream>
#include <armadillo>
#include "new_mesh/ngh/include/polygon.h"
#include "new_mesh/ngh/include/problem.h"

using namespace mathfce;

void GetIntersection(const TBisector & B1, const TBisector &B2,
        TPosition &pos, double &t1, double &t2) {
    TNSolutions ns;
    TMatrix A(2);
    TMVector B(2);
    TMVector X(2);
    double U[ 3 ], V[ 3 ];

    if (AreParallel(B1.GetVector(), B2.GetVector())) {
        if (B1.Belong(B2.GetPoint())) {
            pos = same;
            return;
        } else {
            pos = parallel;
            return;
        }
    }

    B1.GetVector().Get(U);
    B2.GetVector().Get(V);

    int r1 = -1;
    int r2 = -1;
    for (int i = 0; i < 3; i++) {
        if (!IsZero(U[ i ])) {
            r1 = i;
        }
        if (!IsZero(V[ i ])) {
            r2 = i;
        }
    }

    if (r1 == r2) {
        for (int i = 0; i < 3; i++) {
            if ((i != r1) && !IsZero(U[ i ])) {
                r1 = i;
                break;
            }
            if ((i != r2) && !IsZero(V[ i ])) {
                r2 = i;
                break;
            }
        }
    }

    if ((r1 == -1) || (r2 == -1)) {
        pos = skew;
        return;
    }

    int r3 = -1;
    for (int i = 0; i < 3; i++) {
        if ((r1 != i) && (r2 != i)) {
            r3 = i;
            break;
        }
    }

    A.Set(1, 1, U[ r1 ]);
    A.Set(1, 2, -V[ r1 ]);
    A.Set(2, 1, U[ r2 ]);
    A.Set(2, 2, -V[ r2 ]);
    B.Set(1, B2.GetPoint().Get(r1 + 1) - B1.GetPoint().Get(r1 + 1));
    B.Set(2, B2.GetPoint().Get(r2 + 1) - B1.GetPoint().Get(r2 + 1));

    ns = Gauss(A, &X, B);

    if ((ns == no_solution) || (ns == inf_solutions)) {
        pos = skew;
        return;
    }

    if (IsEqual(U[ r3 ] * X.Get(1) - V[ r3 ] * X.Get(2),
            B2.GetPoint().Get(r3 + 1) - B1.GetPoint().Get(r3 + 1))) {
        t1 = X.Get(1);
        t2 = X.Get(2);
        pos = intersecting;
    } else {
        pos = skew;
    }
}

void GetIntersection(const TBisector & B1, const TBisector &B2,
        TPosition &pos, TPoint *P) {
    double t1, t2;
    GetIntersection(B1, B2, pos, t1, t2);
    if (pos != intersecting)
        return;
    *P = B1.GetPoint(t1);
    return;
}

void GetIntersection(const TAbscissa &A1, const TAbscissa &A2,
        TPosition &pos, TPoint *P) {
    double t1, t2;
    if (!QuickIntersectionTest(A1, A2)) {
        pos = skew;
        return;
    }
    GetIntersection(A1, A2, pos, t1, t2);
    if (pos != intersecting)
        return;
    *P = A1.GetPoint(t1);
    return;
}

void GetIntersection(const TAbscissa &A, const TBisector &B,
        TPosition &pos, TPoint *P) {
    double t1, t2;
    GetIntersection(A, B, pos, t1, t2);
    if (pos != intersecting)
        return;
    *P = A.GetPoint(t1);
    return;
}

void GetIntersection(const TBisector &B, const TAbscissa &A,
        TPosition &pos, TPoint *P) {
    GetIntersection(A, B, pos, P);
    return;
}

void GetIntersection(const TAbscissa &A1, const TAbscissa &A2,
        TPosition &pos, double &t1, double &t2) {
    void (*fGetIntersection)(const TBisector &, const TBisector &,
            TPosition &, double &, double &);
    fGetIntersection = GetIntersection;
    fGetIntersection(A1, A2, pos, t1, t2);
    if (pos == intersecting)
        if (t1 > 1 + epsilon || t1 < 0 - epsilon ||
                t2 > 1 + epsilon || t2 < 0 - epsilon)
            pos = skew;
    return;
}

void GetIntersection(const TAbscissa &A, const TBisector &B,
        TPosition &pos, double &t1, double &t2) {
    void (*fGetIntersection)(const TBisector &, const TBisector &,
            TPosition &, double &, double &);
    fGetIntersection = GetIntersection;
    fGetIntersection(A, B, pos, t1, t2);
    if (pos == intersecting)
        if (t1 > 1 + epsilon || t1 < 0 - epsilon)
            pos = skew;
    return;
}

void GetIntersection(const TBisector &B, const TAbscissa &A,
        TPosition &pos, double &t2, double &t1) {
    GetIntersection(A, B, pos, t1, t2);
    return;
}

double Distance(const TBisector & B, const TPoint &P) {
    double d;
    TVector U(P, B.GetPoint());
    d = Cross(U, B.GetVector()).Length() / B.GetVector().Length();
    return d;
}

double Distance(const TPlain &P, const TPoint &X) {
    double dis;
    dis = fabs(P.GetA() * X.X() + P.GetB() * X.Y() + P.GetC() * X.Z() + P.GetD()) /
            sqrt(pow(P.GetA(), 2) + pow(P.GetB(), 2) + pow(P.GetC(), 2));
    return dis;
}

double Distance(const TPoint &P1, const TPoint &P2) {
    double dis;
    dis = sqrt(pow(P1.X() - P2.X(), 2) + pow(P1.Y() - P2.Y(), 2) +
            pow(P1.Z() - P2.Z(), 2));
    return dis;
}

void GetIntersection(const TPlain &P1, const TPlain &P2,
        TPosition &pos, TBisector *B) {
    TVector U;
    TMatrix M(3, 3);
    TMVector b(3);
    TMVector x(3);
    int u1, v1, u2, v2;
    int i, cit;
    bool test;
    if (AreParallel(P1.GetNormal(), P2.GetNormal())) {
        if (P1.Belong(P2.GetPoint())) {
            pos = same;
            return;
        } else {
            pos = parallel;
            return;
        }
    }
    U = Cross(P1.GetNormal(), P2.GetNormal());
    B->SetVector(U);
    for (i = 1; i <= 3; i++)
        b.Set(i, P2.GetPoint().Get(i) - P1.GetPoint().Get(i));
    test = false;
    cit = 0;
    while (!test) {
        u1 = v1 = u2 = v2 = -1;
        switch (cit) {
            case 0: u1 = 1;
                u2 = 2;
                v1 = 3;
                break;
            case 1: u1 = 1;
                u2 = 2;
                v2 = 3;
                break;
            case 2: u2 = 1;
                v1 = 2;
                v2 = 3;
                break;
            case 3: u1 = 1;
                v1 = 2;
                v2 = 3;
                break;
            default:
                if (P1.Belong(P2.GetPoint())) {
                    pos = same;
                    return;
                } else {
                    pos = parallel;
                    return;
                }
                //        mythrow("The two planes should have intersection bisector",__LINE__, __FUNC__);
        }
        for (i = 1; i <= 3; i++) {
            if (u1 != -1)
                M.Set(i, u1, P1.GetU().Get(i));
            if (v1 != -1)
                M.Set(i, v1, P1.GetV().Get(i));
            if (u2 != -1)
                M.Set(i, u2, -P2.GetU().Get(i));
            if (v2 != -1)
                M.Set(i, v2, -P2.GetV().Get(i));
        }
        if (Gauss(M, &x, b) == one_solution)
            test = true;
        cit++;
    }
    if (u1 != -1 && v1 != -1)
        B->SetPoint(P1.GetPoint(x.Get(u1), x.Get(v1)));
    else
        B->SetPoint(P2.GetPoint(x.Get(u2), x.Get(v2)));
    pos = intersecting;
    return;
}

void GetIntersection(const TPlain &P, const TBisector &B,
        TPosition &pos, double &t) {
    arma::mat aa(3,3);
    arma::vec xx(3);
    arma::vec bb(3);
    //TMatrix M(3, 3);
    //TMVector b(3);
    //TMVector x(3);

    int i;
    xx.zeros();
    if (ArePerpendicular(P.GetNormal(), B.GetVector())) {
        if (P.Belong(B.GetPoint())) {
            pos = belong;
            return;
        } else {
            pos = parallel;
            return;
        }
    }
    for (i = 1; i <= 3; i++) {
        /*M.Set(i, 1, B.GetVector().Get(i));
        M.Set(i, 2, -P.GetU().Get(i));
        M.Set(i, 3, -P.GetV().Get(i));
        b.Set(i, P.GetPoint().Get(i) - B.GetPoint().Get(i));*/

        aa(i-1,0) = (B.GetVector().Get(i));
        aa(i-1,1) = (-P.GetU().Get(i));
        aa(i-1,2) = (-P.GetV().Get(i));
        bb(i-1) = (P.GetPoint().Get(i) - B.GetPoint().Get(i));
    }
    arma::solve(xx, aa, bb);
    pos = intersecting;
    t = xx(0);
    return;
}

void GetIntersection(const TPlain &P, const TBisector &B,
        TPosition &pos, TPoint *Pt) {
    double t;
    GetIntersection(P, B, pos, t);
    *Pt = B.GetPoint(t);
}

void GetIntersection(const TBisector &B, const TPlain &P,
        TPosition &pos, TPoint *Pt) {
    GetIntersection(P, B, pos, Pt);
    return;
}

void GetIntersection(const TTriangle &T1, const TTriangle &T2,
        TIntersectionType &it, double &value) {
    //******************************************************************************
    //ONE SHOULD ADD TO THIS FUNCTION POINTS WHICH ARE INSIDE TRIANGLE
    //******************************************************************************
    if (!QuickIntersectionTest(T1, T2)) {
        it = none;
        return;
    }

    TPosition pos;
    TBisector b;
    TPolygon pol;
    double t11, t12, t21, t22, t[ 4 ], t1max, t1min, t2max, t2min;
    int i, j, cit;
    //    TPoint P[6];
    GetIntersection(T1.GetPlain(), T2.GetPlain(), pos, &b);
    if (pos == parallel) {
        it = none;
        return;
    }
    if (pos == intersecting) {
        it = line;
        GetIntersection(b, T1, it, t11, t12);
        if (it == none)
            return;
        GetIntersection(b, T2, it, t21, t22);
        if (it == none)
            return;
        if (t11 > t12) {
            t1max = t11;
            t1min = t12;
        } else {
            t1max = t12;
            t1min = t11;
        }
        if (t21 > t22) {
            t2max = t21;
            t2min = t22;
        } else {
            t2max = t22;
            t2min = t21;
        }
        if (t1max < t2min || t2max < t1min) {
            it = none;
            return;
        }
        t[ 0 ] = t11;
        t[ 1 ] = t12;
        t[ 2 ] = t21;
        t[ 3 ] = t22;
        SortAsc(t, 4);
        value = Distance(b.GetPoint(t[ 1 ]), b.GetPoint(t[ 2 ]));
        return;
    }
    if (pos == same) {
        for (i = 1; i <= 3; i++) {
            cit = 0;
            for (j = 1; j <= 3; j++) {
                GetIntersection(T1.GetAbscissa(i), T2.GetAbscissa(j), pos, t11, t22);
                if (pos == intersecting) {
                    t[ cit ] = t11;
                    cit++;
                }
            }
            if (cit > 1)
                SortAsc(t, 2);
            for (j = 0; j < cit; j++)
                pol.Add(T1.GetAbscissa(i).GetPoint(t[ j ]));
        }
        for (i = 1; i <= 3; i++) {
            if (T1.IsInner(T2.GetPoint(i)))
                pol.Add(T2.GetPoint(i));
            if (T2.IsInner(T1.GetPoint(i)))
                pol.Add(T1.GetPoint(i));
        }
        it = area;
        value = pol.GetArea();
        return;
    }
}

//******************************************************************************
// The value of it is important for computation
//******************************************************************************

void GetIntersection(const TBisector &B, const TTriangle &T,
        TIntersectionType &it, double &t1, double &t2) {

    if (it == none) {
        return;
    }

    TPosition pos;
    double t;
    if (it == unknown || it == point) {
        GetIntersection(T.GetPlain(), B, pos, t);
        switch (pos) {
            case belong: it = line;
                break;
            case intersecting: it = point;
                break;
            default: it = none;
                return;
        }
    }

    if (it == point) {
    	if (!T.IsInner(B.GetPoint(t))) {
            it = none;
            return;
        }
        t1 = t;
        return;
    }

    int cit;
    double tt1, tt2, tt[2];
    if (it == line) {
        cit = 0;
        for (int i = 1; i <= 3; i++) {
            GetIntersection(T.GetAbscissa(i), B, pos, tt1, tt2);
            if (pos == intersecting) {
                tt[ cit ] = tt2;
                cit++;
                if (cit == 2 && IsEqual(tt[ 0 ], tt[ 1 ])) cit = 1;
            }
            if (cit == 2 && !IsEqual(tt[ 0 ], tt[ 1 ]))
                break;
        }
        if (cit != 2) {
            it = none;
            return;
        }
        t1 = tt[ 0 ];
        t2 = tt[ 1 ];
        return;
    }
}

//******************************************************************************
// The value of it is important for computation
//******************************************************************************

void GetIntersection(const TAbscissa &A, const TTriangle &T,
        TIntersectionType &it, double &t1, double &t2) {

	void (*fGetIntersection)(const TBisector &, const TTriangle &,
            TIntersectionType &, double &, double &);

    if (!QuickIntersectionTest(A, T)) {
        it = none;
        return;
    }

    fGetIntersection = GetIntersection;
    fGetIntersection(A, T, it, t1, t2);

    if (it == point || it == line) {
        if (t1 > 1 + epsilon || t1 < 0 - epsilon ||
                t2 > 1 + epsilon || t2 < 0 - epsilon) {
            it = none;
        }
    }

    return;
}

void GetIntersection(const TAbscissa &A, const TTetrahedron &T,
        TIntersectionType &it, double &t1, double &t2) {

    if (!QuickIntersectionTest(A, T)) {
        it = none;
        return;
    }

    int cit = 0;
    double tt1, tt2;
    for (int i = 1; i <= 4; i++) {
        it = unknown;
        GetIntersection(A, T.GetTriangle(i), it, tt1, tt2);
        if (it == line) {
            t1 = tt1;
            t2 = tt2;

            return;
        }
        if (it == point) {
            if (cit == 0) {
                t1 = tt1;
                cit++;
            } else {
                if (IsEqual(t1, tt1)) {
                    continue;
                }
                t2 = tt1;
                it = line;
                return;
            }
        }
    }

    it = none;

    return;
}

void GetIntersection(const TTriangle &Tr, const TTetrahedron &Te,
        TIntersectionType &it, double &coef) {
    //******************************************************************************
    //ONE SHOULD ADD TO THIS FUNCTION POINTS WHICH ARE INSIDE TETRAHEDRON
    //******************************************************************************
    if (!QuickIntersectionTest(Tr, Te)) {
        it = none;
        return;
    }

    TPolygon *P = new TPolygon();

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 4; j++) {
            it = unknown;
            double t1, t2;
            GetIntersection(Tr.GetAbscissa(i), Te.GetTriangle(j), it, t1, t2);
            switch (it) {
                case point:
                    P->Add(Tr.GetAbscissa(i).GetPoint(t1));
                    break;
                case line:
                    P->Add(Tr.GetAbscissa(i).GetPoint(t1));
                    P->Add(Tr.GetAbscissa(i).GetPoint(t2));
                    break;
                default:
                    //mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
                	break;
            }
        }
    }

    for (int i = 1; i <= 6; i++) {
        it = unknown;
        double t1, t2;
        GetIntersection(Te.GetAbscissa(i), Tr, it, t1, t2);
        switch (it) {
            case point:
                P->Add(Te.GetAbscissa(i).GetPoint(t1));
                break;
            case line:
                P->Add(Te.GetAbscissa(i).GetPoint(t1));
                P->Add(Te.GetAbscissa(i).GetPoint(t2));
                break;
            default:
                //mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
            	break;
        }
    }

    coef = P->GetArea();
    if (IsZero(coef)) {
        it = none;
    } else {
        it = area;
    }

    return;
}

//******************************************************************************
// This function returns true if two input geometry can intersect
//******************************************************************************

template<class A, class B> bool QuickIntersectionTest(const A &a, const B &b) {
    for (int i = 1; i <= 3; i++) {
        if (a.GetMin(i) > b.GetMax(i) || a.GetMax(i) < b.GetMin(i)) {
            return false;
        }
    }
    return true;
}



