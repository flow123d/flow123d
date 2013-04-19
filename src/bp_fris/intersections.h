#include "simplex.h"

// Vytvořit objekty typu Intersections


SPoint<3> Intersections(HyperPlane<1,3> a, Simplex<2,3> b);

SPoint<3> Intersections(HyperPlane<1,3> a, Simplex<1,3> b);

void Intersections(HyperPlane<1,3> a, Simplex<3,3> b, SPoint<3> &prus1, SPoint<3> &prus2);

bool Intersections(HyperPlane<1,3> a, Simplex<0,3> b);

SPoint<3>* IntersectionsOpp(HyperPlane<1,3> hp, Simplex<2,3> sm, int stena);

void IntersectionsOp(HyperPlane<1,3> hp, Simplex<3,3> sm, SPoint<3> &prusecik1, SPoint<3> &prusecik2);

SPoint<3> Globalni_souradnice(SPoint<3> lokalni_souradnice, Simplex<2,3> trojuhelnik);

bool IntersectionsOp_1D_2D(HyperPlane<1,3> &hp, Simplex<2,3> sm, int stena, std::vector<double> &coords_3D, double &local_abscissa, bool &orientace);
