#ifndef intersectionH
#define intersectionH

#include "bisector.h"
#include "abscissa.h"
#include "point.h"
#include "plain.h"
#include "triangle.h"
#include "tetrahedron.h"
#include "intersectionLocal.h"

typedef enum Intersections{
        none,
        unknown,
        point,
        line,
        area
} TIntersectionType;

typedef enum Positions{
        skew,
        parallel,
        intersecting,
        same,
        belong
} TPosition;

void GetIntersection(const TBisector &, const TBisector &, TPosition &,
                     double &, double &);
//void GetIntersection(const TBisector &, const TBisector &, TPosition &,
//                     TPoint *);
void GetIntersection(const TAbscissa &, const TAbscissa &, TPosition &,
                     double &, double &);
//void GetIntersection(const TAbscissa &, const TAbscissa &, TPosition &, //puvodni
//                     TPoint *);
void GetIntersection(const TAbscissa &, const TAbscissa &, IntersectionLocal * & insec);

void GetIntersection(const TBisector &, const TAbscissa &, TPosition &,
                     double &, double &);
//void GetIntersection(const TBisector &, const TAbscissa &, TPosition &, //puvodni
//                     TPoint *);
void GetIntersection(const TBisector &, const TAbscissa &, IntersectionLocal * & insec);

void GetIntersection(const TAbscissa &, const TBisector &, TPosition &,
                     double &, double &);
//void GetIntersection(const TAbscissa &, const TBisector &, TPosition &, //puvodni
//                     TPoint *);
void GetIntersection(const TAbscissa &, const TBisector &, IntersectionLocal * & insec);

void GetIntersection(const TPlain &, const TPlain &,
                     TPosition &, TBisector *);
void GetIntersection(const TPlain &, const TBisector &,
                     TPosition &, TPoint *);
void GetIntersection(const TBisector &, const TPlain &,
                     TPosition &, double &);
void GetIntersection(const TBisector &, const TPlain &,
                     TPosition &, TPoint *);
void GetIntersection(const TTriangle &, const TTriangle &,
                     TIntersectionType &, double &);
//void GetIntersection(const TBisector &, const TTriangle &, //puvodni
//                     TIntersectionType &, double &, double &);
void GetIntersection(const TBisector &, const TTriangle &, IntersectionLocal * & insec);

//void GetIntersection(const TAbscissa &, const TTriangle &, //puvodni
//                     TIntersectionType &, double &, double &);
void GetIntersection(const TAbscissa &, const TTriangle &, IntersectionLocal * & insec);

void GetIntersection(const TAbscissa &, const TTetrahedron &,
                     TIntersectionType &, double &);
void GetIntersection(const TTriangle &, const TTetrahedron &,
                     TIntersectionType &, double &);

template<class A, class B> bool QuickIntersectionTest(const A &a, const B &b);

double Distance(const TBisector &, const TPoint &);
double Distance(const TPlain &, const TPoint &);
double Distance(const TPoint &, const TPoint &);
#endif
