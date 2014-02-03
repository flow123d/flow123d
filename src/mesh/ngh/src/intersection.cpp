#include "system/exc_common.hh"
#include "mesh/ngh/include/intersection.h"
#include "mesh/ngh/include/intersectionLocal.h"
#include "mesh/ngh/include/matrix.h"
#include "mesh/ngh/include/mathfce.h"
#include <cmath>
//#include "math.h"
#include <iostream>
#include <armadillo>
#include "mesh/ngh/include/polygon.h"
//#include "mesh/ngh/include/problem.h"

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
/*
void GetIntersection(const TBisector & B1, const TBisector &B2,
        TPosition &pos, TPoint *P) {
    double t1, t2;
    GetIntersection(B1, B2, pos, t1, t2);
    if (pos != intersecting)
        return;
    *P = B1.GetPoint(t1);
    return;
}
*/
/*
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
*/

//------------------------UPRAVENO---------------------------
// vraci localni souradnice pruniku, prvni souradnice vzhledem k A1 druha souradnice vzhledem k A2

void GetIntersection(const TAbscissa &A1, const TAbscissa &A2, IntersectionLocal * &insec) {
	double t1, t2;
	TPosition pos;

    if (!QuickIntersectionTest(A1, A2)) {
        insec=NULL;
        return;
    }
    GetIntersection(A1, A2, pos, t1, t2);
    if ( pos == intersecting ) {
        // test t1 je (0-eps,1+eps) a t2 je z (0-eps,1+eps)
    	if ((t1 > (0 - epsilon)) || (t1 < (1 + epsilon)) || (t2 > (0 - epsilon)) || (t2 < (1 + epsilon))) {
    		insec=new IntersectionLocal(IntersectionLocal::point);
    		vector<double> loc_coord_1(1,t1);
    		vector<double> loc_coord_2(1,t2);
    		insec->add_local_coord(loc_coord_1,loc_coord_2);
    		return;
    	} else {
    		insec = NULL;
    		return;
    	}
    } else if ( pos == same ) { //A1 a A2 lezi na stejne primce
        // nalezeni prunikove usecky X1, X2 ; Y1, Y2  ...koncove body usecky XY

    	TAbscissa A2_copy(A2);
    	double dot_product;
    	// TEST - vzajemna orientace A1 a A2 (skalarni soucin vektoru A1, A2)
    	dot_product = Dot(A1.GetVector(), A2.GetVector());
    		if (dot_product < 0 ) { //if opacna orientace A1 a A2
    			TVector AX2((-1)*A2.GetVector());
    			A2_copy.SetVector(AX2);
    		}

    //vypocet lokalni souradnice koncu A2 vzhledem k A1
    	TVector Diff;
    	double loc_begin_A1, loc_end_A1, loc_begin_A2, loc_end_A2;
    	Diff=A2_copy.GetPoint()-A1.GetPoint();
    	loc_begin_A2=Diff.Length()/(A1.GetVector().Length());

    	Diff=A2_copy.GetPoint(1)-A1.GetPoint();
    	loc_end_A2=Diff.Length()/(A1.GetVector().Length());

    //vypocet lokalni souradnice koncu A1 vzhledem k A2
    	Diff=A1.GetPoint()-A2_copy.GetPoint();
    	loc_begin_A1=Diff.Length()/(A2_copy.GetVector().Length());

    	Diff=A1.GetPoint(1)-A2_copy.GetPoint();
    	loc_end_A1=Diff.Length()/(A2_copy.GetVector().Length());

    //X1...loc.souradnice X vzhledem k A1
    //X2...loc.souradnice X vzhledem k A2
    //Y1...loc.souradnice Y vzhledem k A1
    //Y2...loc.souradnice Y vzhledem k A2
    	double X1, X2, Y1, Y2;
    //1.possibility - A2 lezi pred A1
    	if ((loc_begin_A2 < 0) && (loc_end_A2 < 0)) {
    		insec=NULL;
    		return;
    	}
   	//2.possibility - prekryvaji se, A2 lezi pred A1
    	if ((loc_begin_A2 < 0) && (loc_end_A2 > 0) && (loc_end_A2 < 1)) {
    		X1 = 0;
    		X2 = loc_begin_A1;
    		Y1 = loc_end_A2;
    		Y2 = 1;
    	}
    //3.possibility - prekryvaji se, A2 lezi mezi koncovymi body A1
    	if ((loc_begin_A2 > 0) && (loc_end_A2 > 0) && (loc_end_A2 < 1)) {
    		X1 = loc_begin_A2;
    		X2 = 0;
    		Y1 = loc_end_A2;
    		Y2 = 1;
    	}
    //4.possibility - prekryvaji se, A2 lezi za A1
    	if ((loc_begin_A2 > 0) && (loc_begin_A2 < 1) && (loc_end_A2 > 1)) {
    		X1 = loc_begin_A2;
    		X2 = 0;
    		Y1 = 1;
    		Y2 = loc_end_A1;
    	}
    //5.possibility - A2 lezi za A1
    	if ((loc_begin_A2 > 1) && (loc_end_A2 > 1)) {
    		insec=NULL;
    		return;
    	}
    //6.possibility - prekryvaji se, A1 lezi mezi koncovymi body A2
    	if ((loc_begin_A2 < 0) && (loc_end_A2 > 1)) {
    		X1 = 0;
    		X2 = loc_begin_A1;
    		Y1 = 1;
    		Y2 = loc_end_A1;
    	}

    //set local coords:
    	insec=new IntersectionLocal(IntersectionLocal::line);
    	vector<double> loc_coord_1(1);
    	vector<double> loc_coord_2(1);
    	loc_coord_1[0] = X1;
    	loc_coord_2[0] = X2;
    	insec->add_local_coord(loc_coord_1,loc_coord_2);
    	loc_coord_1[0] = Y1;
    	loc_coord_2[0] = Y2;
    	insec->add_local_coord(loc_coord_1,loc_coord_2);
    	return;

    } else {
        insec = NULL;
        return;
    }
    insec = NULL;
	return;
}
/*
void GetIntersection(const TAbscissa &A, const TBisector &B,
        TPosition &pos, TPoint *P) {
    double t1, t2;
    GetIntersection(A, B, pos, t1, t2);
    if (pos != intersecting)
        return;
    *P = A.GetPoint(t1);
    return;
}
*/

//------------------------UPRAVENO---------------------------
// vraci localni souradnice pruniku, prvni souradnice vzhledem k A druha souradnice vzhledem k B
// v pripade ze je prunik usecka - poradi bodu podle orientace abscissa A

void GetIntersection(const TAbscissa &A, const TBisector &B, IntersectionLocal * &insec) {
    double t1, t2;
    TPosition pos;
    insec = NULL;

    GetIntersection(A, B, pos, t1, t2);
    if ( pos == intersecting ) {
        // test t1 je (0-eps,1+eps)
        if ((t1 > (0 - epsilon)) || (t1 < (1 + epsilon))) {
        	insec=new IntersectionLocal(IntersectionLocal::point);
        	vector<double> loc_coord_1(1,t1);
        	vector<double> loc_coord_2(1,t2); //t2 na Bisectoru B
        	insec->add_local_coord(loc_coord_1,loc_coord_2);
        }
    } else if ( pos == same ) { //A a B lezi na stejne primce
        // prunik cela usecka =>ulozit koncove body Abscissy A

    //vypocet lokalni souradnice koncovych bodu A
    	bool x;

    	// get local coords of points of A on B
    	B.GetParameter(A.GetPoint(), t1, x);
    	B.GetParameter(A.GetPoint(1), t2, x);

    	// sort t1, t2
    	if (t2 < t1) {
    		double swap = t1;
    		t1 = t2;
    		t2 = swap;
    	}

    	//set local coords:
    	if ((t2 > (0 - epsilon)) && (t1 < (1 + epsilon))) {
    		double loc_begin_A, loc_end_A;
    		loc_begin_A = (t1 > 0) ? t1 : 0;
    		loc_end_A = (t2 < 1) ? t2 : 1;

    		insec=new IntersectionLocal(IntersectionLocal::line);
			// first abscissa point
			vector<double> loc_coord_1(1,0);
			vector<double> loc_coord_2(1,loc_begin_A);
			insec->add_local_coord(loc_coord_1,loc_coord_2);
			// second abscissa point
			vector<double> loc_coord_1_(1,1);
			vector<double> loc_coord_2_(1,loc_end_A);
			insec->add_local_coord(loc_coord_1_,loc_coord_2_);
			return;
    	}
    } else {
    	insec=NULL;
    	return;
    }
    return;
}

/*
void GetIntersection(const TBisector &B, const TAbscissa &A,
        TPosition &pos, TPoint *P) {
    GetIntersection(A, B, pos, P);
    return;
}
*/
void GetIntersection(const TBisector &B, const TAbscissa &A, IntersectionLocal * &insec) {
    GetIntersection(A, B, insec); //KONTROLA
    return;
}

void GetIntersection(const TAbscissa &A1, const TAbscissa &A2,
        TPosition &pos, double &t1, double &t2) {

    GetIntersection( (const TBisector &)A1, (const TBisector &)A2, pos, t1, t2);
    if (pos == intersecting)
        if (t1 > 1 + epsilon || t1 < 0 - epsilon ||
                t2 > 1 + epsilon || t2 < 0 - epsilon)
            pos = skew;
    return;
}

void GetIntersection(const TAbscissa &A, const TBisector &B,
        TPosition &pos, double &t1, double &t2) {

    GetIntersection( (const TBisector &)A, B, pos, t1, t2);
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
        TIntersectionType &it, double &value) {	//ZATIM NEPOTREBUJEME =>zakomentovano

    //******************************************************************************
    //ONE SHOULD ADD TO THIS FUNCTION POINTS WHICH ARE INSIDE TRIANGLE
    //******************************************************************************
/*
    if (!QuickIntersectionTest(T1, T2)) {
        it = none;
        return;
    }

    TPosition pos;
    TPoint P1,P2;
    TBisector b(P1,P2);
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
*/
}

//******************************************************************************
// The value of it is important for computation
//******************************************************************************

//------------------------UPRAVENO---------------------------
// vraci localni souradnice pruniku, prvni souradnice vzhledem k B druha souradnice vzhledem k T
// pro vsechny body pruniku B a T

void GetIntersection(const TBisector &B, const TTriangle &T, IntersectionLocal * &insec) {
    TPosition pos;
    double t;
    int cit=0;
    GetIntersection(T.GetPlain(), B, pos, t);
    switch (pos) {
        // POINT INTERSECTION
    	case intersecting: {
            if (!T.IsInner(B.GetPoint(t))) {
                insec=NULL;
                return;
            }

            vector<double> loc_coord_1(1);
            vector<double> loc_coord_2(2);

            loc_coord_1[0] = t;
            loc_coord_2[0] = 0; //TODO: values of loc_coord_2 are not computed
            loc_coord_2[1] = 0;

            insec = new IntersectionLocal(IntersectionLocal::point);
            insec->add_local_coord(loc_coord_1, loc_coord_2);

            return;
    	}

       // LINE INTERSECTION
       case belong: {
           IntersectionLocal* insec_tmp;
           IntersectionPoint* insec_point_tmp[3];

           vector<double> loc_tria_tmp(2);
           // inicializace vektoru pro lokalni souradnice pro pripad useckoveho pruseciku
           // lokalni souradnice vrholu trojuhelnika
           vector<double> loc_tria_coord_01(2); //loc. coord 1.vrcholu trojuhelnika
           vector<double> loc_tria_coord_02(2); //loc. coord 2.vrcholu trojuhelnika
           vector<double> loc_tria_coord_03(2); //loc. coord 3.vrcholu trojuhelnika

           loc_tria_coord_01[0]=0;
           loc_tria_coord_01[1]=0;

           loc_tria_coord_02[0]=1;
           loc_tria_coord_02[1]=0;

           loc_tria_coord_03[0]=0;
           loc_tria_coord_03[1]=1;

           //PRVNI STENA TROJUHELNIKU
           GetIntersection(T.GetAbscissa(1), B, insec_tmp);
           if (insec_tmp != NULL) {
        	   if (insec_tmp->get_type() == IntersectionLocal::point) {
        		   // intersection in a point
        		   // zde: 1) vybrat IntersectionPoint z insec_tmp
                   //      2) prohodit jeho coord
        		   //      3) z druhe coord udelat lokalni souradnici na trojuhelniku

        		   // loc. triangle coord (X, 0)
                   loc_tria_tmp[0] = insec_tmp->get_point(0)->el1_coord()[0];
                   loc_tria_tmp[1] = 0;
        		   insec_point_tmp[cit] = new IntersectionPoint(insec_tmp->get_point(0)->el2_coord(), loc_tria_tmp);
        		   cit++; //citac kolik sten protne
        		   //cout<<"\nCit(1.stena)= "<< cit << endl;
        		   delete insec_tmp;
               } else if (insec_tmp->get_type() == IntersectionLocal::line) {
            	   // intersection is whole tringle side => ulozit lokalni souradnice strany trojuhelnika do insec
                   insec = new IntersectionLocal(IntersectionLocal::line);
                   insec->add_local_coord(insec_tmp->get_point(0)->el2_coord(), loc_tria_coord_01);
                   insec->add_local_coord(insec_tmp->get_point(1)->el2_coord(), loc_tria_coord_02);
                   delete insec_tmp;
                   return;
               }
           }

           //DRUHA STENA TROJUHELNIKU
           GetIntersection(T.GetAbscissa(2), B, insec_tmp);
           if (insec_tmp != NULL) {
        	   if (insec_tmp->get_type() == IntersectionLocal::point) {
        		   // loc. triangle coord (1-X, X)
        		   loc_tria_tmp[0] = 1 - insec_tmp->get_point(0)->el1_coord()[0];
        		   loc_tria_tmp[1] = insec_tmp->get_point(0)->el1_coord()[0];
        		   insec_point_tmp[cit] = new IntersectionPoint(insec_tmp->get_point(0)->el2_coord(), loc_tria_tmp);
        		   cit++; //citac kolik sten protne
        		   //cout<<"Cit(2.stena)= "<< cit << endl;
         		   delete insec_tmp;
               } else if (insec_tmp->get_type() == IntersectionLocal::line) {
            	   // intersection is whole tringle side => ulozit lokalni souradnice strany trojuhelnika do insec
                   insec = new IntersectionLocal(IntersectionLocal::line);
                   insec->add_local_coord(insec_tmp->get_point(0)->el2_coord(), loc_tria_coord_02);
                   insec->add_local_coord(insec_tmp->get_point(1)->el2_coord(), loc_tria_coord_03);
                   delete insec_tmp;
                   return;
               }
           }

           //TRETI STENA TROJUHELNIKU
           GetIntersection(T.GetAbscissa(3), B, insec_tmp);
           if (insec_tmp != NULL) {
        	   if (insec_tmp->get_type() == IntersectionLocal::point) {
        		   // loc. triangle coord (0, 1-X)
        		   loc_tria_tmp[0] = 0;
        		   loc_tria_tmp[1] = 1 - insec_tmp->get_point(0)->el1_coord()[0];
        		   insec_point_tmp[cit] = new IntersectionPoint(insec_tmp->get_point(0)->el2_coord(), loc_tria_tmp);
        		   cit++; //citac kolik sten protne
        		   //cout <<"Cit(3.stena)= "<< cit << endl;
        		   delete insec_tmp;
               } else if (insec_tmp->get_type() == IntersectionLocal::line) {
            	   // intersection is whole tringle side => ulozit lokalni souradnice strany trojuhelnika do insec
                   insec = new IntersectionLocal(IntersectionLocal::line);
                   insec->add_local_coord(insec_tmp->get_point(0)->el2_coord(), loc_tria_coord_03);
                   insec->add_local_coord(insec_tmp->get_point(1)->el2_coord(), loc_tria_coord_01);
                   delete insec_tmp;
                   return;
               }
           }

           //TEST - pocet protnutych sten trojuhelnika
           if (cit == 0) {
               insec = NULL;
               return;
           }
           //lezi pres vrchol a zaroven protne protilehlou stenu => cit==3
           if (cit == 3) {
        	   //cout<<"(cit == 3) => REDUKCE insec_point_tmp!" << endl;
        	   if (*(insec_point_tmp[0]) == *(insec_point_tmp[0])) {
        		   //cout<<"2 insec_point_tmp[0] jsou stejne" << endl;
        	   }

        	   for (int i = 0; i < cit; i++) {
        		   if (*(insec_point_tmp[i]) == *(insec_point_tmp[(i+1)%3])) {

        			   //cout<<"insec_point_tmp[(i+1)%3] :: el1_coord().size(): " << insec_point_tmp[(i+1)%3]->el1_coord().size() << endl;
        			   //cout<<"(insec_point_tmp[(i+1)%3]) :: el2_coord().size(): " << (insec_point_tmp[(i+1)%3])->el2_coord().size() << endl;

        			   delete insec_point_tmp[(i+1)%3];
        			   cit--;
        			   //cout<<"PO REDUKCI cit= " << cit << endl;

        			   //nutno posunout zbyle prvky pole insec_point_tmp
        			   IntersectionPoint* tmp = insec_point_tmp[(i+2)%3];
        			   insec_point_tmp[0] = insec_point_tmp[i];
        			   insec_point_tmp[1] = tmp;
        			   break;
        		   }
        	   }
           }
           if (cit != 2) {
        	   THROW( ExcAssertMsg() << EI_Message("Error - pocet bodu pruniku != 2.\n") );
        	   return;
           } else {
        	   if (*(insec_point_tmp[0]) == *(insec_point_tmp[1])) { //lezi pres vrchol
        		   insec = new IntersectionLocal(IntersectionLocal::point);
        		   insec->add_local_point(insec_point_tmp[0]);
        		   //cout<<"Insec_point_1 == Insec_point_2" << endl;
        		   delete insec_point_tmp[1];
        	   } else {
        		   insec = new IntersectionLocal(IntersectionLocal::line);
        		   insec->add_local_point(insec_point_tmp[0]);
        		   insec->add_local_point(insec_point_tmp[1]);
        		   //cout<<"(IntersectionType=line) - point_1->el1_coord().size()= " << insec->get_point(0)->el1_coord().size()<< endl;
        		   //cout<<"(IntersectionType=line) - point_2->el1_coord().size()= " << insec->get_point(1)->el1_coord().size()<< endl;
        	   }
        	   return;
           }
        } //end case belong

        // EMPTY INTERSECTION
        default:
        	insec = NULL;
        	return;
        }
}

//******************************************************************************
// The value of it is important for computation
//******************************************************************************

//------------------------UPRAVENO---------------------------
// vraci localni souradnice pruniku, prvni souradnice vzhledem k A druha souradnice vzhledem k T
// pro vsechny body pruniku A a T

void GetIntersection(const TAbscissa &A, const TTriangle &T,
        IntersectionLocal * & insec) {

    if (!QuickIntersectionTest(A, T)) {
        insec = NULL;
        return;
    }
    IntersectionLocal* insec_tmp;
    GetIntersection( (const TBisector &)A, T, insec_tmp);
    if (!insec_tmp) {
    	insec = NULL;
    	return;
    }
    if (insec_tmp->get_type() == IntersectionLocal::point) {
    	if (insec_tmp->get_point(0) != NULL) {
    		double t1 = insec_tmp->get_point(0)->el1_coord()[0];
    		if (t1 < 0 - epsilon || t1 > 1 + epsilon) {
    		    delete insec_tmp;
    			insec = NULL;
    		} else {
    		    insec = insec_tmp;
    		}
    	}
    } else if(insec_tmp->get_type() == IntersectionLocal::line) {
        // A1 i A2 ma byt v intervalu (0,1) -> vrati insec
        // pokud ne tak zkusi zkratit, nebo NULL (delete)

    	IntersectionPoint* A1;
    	IntersectionPoint* A2;
    	if (insec_tmp->get_point(0) != NULL) {
    	    if (insec_tmp->get_point(0)->el1_coord()[0] > insec_tmp->get_point(1)->el1_coord()[0]) {
    	    	A2 = insec_tmp->get_point(0);
    	    	A1 = insec_tmp->get_point(1);
    	    } else {
    	    	A1 = insec_tmp->get_point(0);
    	    	A2 = insec_tmp->get_point(1);
    	    }

    	    double A1_t = A1->el1_coord()[0];
    	    double A2_t = A2->el1_coord()[0];
    	    if (A1_t < 0) A1_t = 0;
    	    if (A2_t > 1) A2_t = 1;

    	    if (A2_t < A1_t) {
    			delete insec_tmp;
    			insec = NULL;
    	    } else {
				insec = new IntersectionLocal(IntersectionLocal::line);
				//cout << "A1_t: " << A1_t << endl;
				//cout << "A2_t: " << A2_t << endl;
				insec->add_local_point(interpolate(*A1, *A2, A1_t));
				insec->add_local_point(interpolate(*A1, *A2, A2_t));
				delete insec_tmp;
    	    }
    	}
    return;
    }
    return;

//-------------------------------------------------------------------
/*    		IntersectionPoint* insec_point_tmp[2];
    		vector<double> loc_tria_tmp(2);
    		double A1_ = insec_tmp->get_point(0)->el1_coord()[0];
    		double A2_ = insec_tmp->get_point(1)->el1_coord()[0]; //test != NULL (lezi pres vrchol!)
    		double *A1 = &A1_;
    		double *A2 = &A2_;
    		bool invert = false;
    		//TEST VZAJEMNE ORIENTACE A1, A2:
    		if (A1_ > A2_) {
    			A1 = &A2_;
    			A2 = &A1_;
    			invert = true;
    		}
    		// pripady:
    		// 1) prusecik uvnitr usecky -> predat nezmenene
    		// 2) prusecik mimo usecku -> delete a vratit NULL
    		// 3) else a) A1 < 0 - zkratit zleva
    		//         b) A2 > 1  - zkratit zprava

    		//PRUSECIK UVNITR USECKY
    		if ((*A1 > 0 - epsilon) && (*A1 < 1 + epsilon) && (*A2 > 0 - epsilon) && (*A2 < 1 + epsilon)) {
    			insec = insec_tmp;
    		}
    		//PRUSECIK MIMO USECKU
    		if (((*A1 < 0 - epsilon) && (*A2 < 0 - epsilon)) || ((*A1 > 1 + epsilon) && (*A2 > 1 + epsilon))) {
    			delete insec_tmp;
    			insec = NULL;
    		} else {
    			//1.possibility (A1 < 0) - zkratit zleva
    			if ((*A1 < 0 - epsilon) && (*A2 > 0 - epsilon) && (*A2 < 1 + epsilon)) {
    				if(invert == true) {
    					//loc. coord. r:
    					loc_tria_tmp[0] = ((1 - *A2)/(*A1 - *A2))*(insec_tmp->get_point(1)->el2_coord()[0] - insec_tmp->get_point(0)->el2_coord()[0]) + insec_tmp->get_point(0)->el2_coord()[0];
    					//loc. coord. s:
    					loc_tria_tmp[1] = ((1 - *A2)/(*A1 - *A2))*(insec_tmp->get_point(1)->el2_coord()[1] - insec_tmp->get_point(0)->el2_coord()[1]) + insec_tmp->get_point(0)->el2_coord()[1];
    				} else {
    					//loc. coord. r:
    					loc_tria_tmp[0] = ((1 - *A2)/(*A1 - *A2))*(insec_tmp->get_point(0)->el2_coord()[0] - insec_tmp->get_point(1)->el2_coord()[0]) + insec_tmp->get_point(1)->el2_coord()[0];
    					//loc. coord. s:
    					loc_tria_tmp[1] = ((1 - *A2)/(*A1 - *A2))*(insec_tmp->get_point(0)->el2_coord()[1] - insec_tmp->get_point(1)->el2_coord()[1]) + insec_tmp->get_point(1)->el2_coord()[1];
    				}
    				vector<double> loc_A1(1, 0);
    				vector<double> loc_A2(1, *A2);
    				if(invert == true) {
    					insec_point_tmp[0] = new IntersectionPoint(loc_A2, insec_tmp->get_point(1)->el2_coord());
    					insec_point_tmp[1] = new IntersectionPoint(loc_A1, loc_tria_tmp);
    				} else {
    					insec_point_tmp[0] = new IntersectionPoint(loc_A2, insec_tmp->get_point(0)->el2_coord());
    					insec_point_tmp[1] = new IntersectionPoint(loc_A1, loc_tria_tmp);
    				}
    				insec = new IntersectionLocal(IntersectionLocal::line);
    				insec->add_local_point(insec_point_tmp[0]);
    				insec->add_local_point(insec_point_tmp[1]);
    				delete insec_tmp;
    			}
    			//2.possibility (A2 > 1) - zkratit zprava
    			if ((*A1 > 0 - epsilon) && (*A1 < 1 + epsilon) && (*A2 > 1 + epsilon)) {
    				if(invert == true) {
    					//loc. coord. r:
    				    loc_tria_tmp[0] = ((1 - *A1)/(*A2 - *A1))*(insec_tmp->get_point(0)->el2_coord()[0] - insec_tmp->get_point(1)->el2_coord()[0]) + insec_tmp->get_point(1)->el2_coord()[0];
    				    //loc. coord. s:
    				    loc_tria_tmp[1] = ((1 - *A1)/(*A2 - *A1))*(insec_tmp->get_point(0)->el2_coord()[1] - insec_tmp->get_point(1)->el2_coord()[1]) + insec_tmp->get_point(1)->el2_coord()[1];
    				} else {
    				    //loc. coord. r:
    				    loc_tria_tmp[0] = ((1 - *A1)/(*A2 - *A1))*(insec_tmp->get_point(1)->el2_coord()[0] - insec_tmp->get_point(0)->el2_coord()[0]) + insec_tmp->get_point(0)->el2_coord()[0];
    				    //loc. coord. s:
    				    loc_tria_tmp[1] = ((1 - *A1)/(*A2 - *A1))*(insec_tmp->get_point(1)->el2_coord()[1] - insec_tmp->get_point(0)->el2_coord()[1]) + insec_tmp->get_point(0)->el2_coord()[1];
    				}
    				vector<double> loc_A1(1, *A1);
    				vector<double> loc_A2(1, 1);
    				if(invert == true) {
    					insec_point_tmp[0] = new IntersectionPoint(loc_A1, insec_tmp->get_point(1)->el2_coord());
    					insec_point_tmp[1] = new IntersectionPoint(loc_A2, loc_tria_tmp);
    				} else {
    					insec_point_tmp[0] = new IntersectionPoint(loc_A1, insec_tmp->get_point(0)->el2_coord());
    					insec_point_tmp[1] = new IntersectionPoint(loc_A2, loc_tria_tmp);
    				}
    				insec = new IntersectionLocal(IntersectionLocal::line);
    				insec->add_local_point(insec_point_tmp[0]);
    				insec->add_local_point(insec_point_tmp[1]);
    				delete insec_tmp;
    			}

    			//3.possibility - zkratit z obou stran
    			//if ((*A1 < 0 - epsilon) && (*A2 > 1 + epsilon)) {
    			//	vector<double> loc_A1(1, 0);
    			//	vector<double> loc_A2(1, 1);
    			//}
    		}
    	}
    return;
    }
    return;
*/
}

void GetIntersection(const TAbscissa &A, const TTetrahedron &T,
        TIntersectionType &it, double &coef) {

    if (!QuickIntersectionTest(A, T)) {
        it = none;
        return;
    }

    int cit = 0;
    double tt1, tt2;
    double tt[2];

    if (T.IsInner(A.GetPoint())) {
    	tt[cit] = 0;
    	cit++;
    }
    if (T.IsInner(A.GetPoint(1))) {
    	tt[cit] = 1;
    	cit++;
    }
    if (cit == 2) {
    	coef = fabs(tt[1] - tt[0]);
    	it = line;
    	return;
    }

    IntersectionLocal *insec;

    for (int i = 1; i <= 4; i++) {
        it = unknown;
        GetIntersection(A, T.GetTriangle(i), insec);
        if (insec) {
            if (insec->get_type() == IntersectionLocal::line) {
            	tt1 = insec->get_point(0)->el1_coord()[0];
            	tt2 = insec->get_point(1)->el1_coord()[0];
            	if (tt1 > tt2) {
            		double swap = tt1;
            		tt1 = tt2;
            		tt2 = swap;
            	}

            	if ((tt2 > (0 - epsilon)) && (tt1 < (1 + epsilon))) {
            		if (tt1 < 0) tt1 = 0;
            		if (tt2 > 1) tt2 = 1;
            		coef = fabs(tt2 - tt1);
            		it = line;
            		return;
            	}
            }
            if (insec->get_type() == IntersectionLocal::point) {
                if (cit == 0) {
                    tt[0] = insec->get_point(0)->el1_coord()[0];
                    cit++;
                } else {
                    if (IsEqual(tt[0], insec->get_point(0)->el1_coord()[0])) {
                        continue;
                    }
                    if (tt[0] > insec->get_point(0)->el1_coord()[0]) {
                    	tt[1] = tt[0];
                    	tt[0] = insec->get_point(0)->el1_coord()[0];
                    } else {
                    	tt[1] = insec->get_point(0)->el1_coord()[0];
                    }
                    if ((tt[1] > (0 - epsilon)) && (tt[0] < (1 + epsilon))) {
						if (tt[0] < 0) tt[0] = 0;
						if (tt[1] > 1) tt[1] = 1;
						coef = fabs(tt[1] - tt[0]);
						it = line;
						return;
                    }
                }
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

    TPolygon P;// = new TPolygon();
    for (int i = 1; i <= 3; i++) {
        if (Te.IsInner(Tr.GetPoint(i))) {
        	P.Add(Tr.GetPoint(i));
        }
    }

    if (P.vertexes_count() < 3) {
		IntersectionLocal *insec;

		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 4; j++) {
				GetIntersection(Tr.GetAbscissa(i), Te.GetTriangle(j), insec);
				if (insec) {
					switch (insec->get_type()) {
						case IntersectionLocal::point: {
							P.Add(Tr.GetAbscissa(i).GetPoint(insec->get_point(0)->el1_coord()[0]));
							break;
						}
						case IntersectionLocal::line: {
							P.Add(Tr.GetAbscissa(i).GetPoint(insec->get_point(0)->el1_coord()[0]));
							P.Add(Tr.GetAbscissa(i).GetPoint(insec->get_point(1)->el1_coord()[0]));
							break;
						}
						default:
							//mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
							break;
					}
				}
			}
		}
		for (int i = 1; i <= 6; i++) {
			GetIntersection(Te.GetAbscissa(i), Tr, insec);
			if (insec) {
				switch (insec->get_type()) {
					case IntersectionLocal::point:
						P.Add(Te.GetAbscissa(i).GetPoint(insec->get_point(0)->el1_coord()[0]));
						break;
					case IntersectionLocal::line:
						P.Add(Te.GetAbscissa(i).GetPoint(insec->get_point(0)->el1_coord()[0]));
						P.Add(Te.GetAbscissa(i).GetPoint(insec->get_point(1)->el1_coord()[0]));
						break;
					default:
						//mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
						break;
				}
			}
		}
    }

    coef = P.GetArea();
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
    	//cout << "i: " << i << "; a: " << a.GetMin(i) << "  "<< a.GetMax(i)<< "; b: " << b.GetMin(i) << "  "<< b.GetMax(i) << endl;
        if (a.GetMin(i) > b.GetMax(i) + epsilon || a.GetMax(i) < b.GetMin(i) - epsilon) {
            return false;
        }
    }
    return true;
}

