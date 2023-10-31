#ifndef FEM_TOOLS_TEST_HH_
#define FEM_TOOLS_TEST_HH_

#include <vector>
#include "system/armor.hh"


static const uint SIZE = 16;

inline double* sum(double* a, double* b, double* c) {
	double *r = new double[SIZE];
    for (uint i=0; i<SIZE; ++i)
        r[i] = a[i] + b[i] + c[i];
    return r;
}

inline double* diff(double* a, double* b) {
	double *r = new double[SIZE];
    for (uint i=0; i<SIZE; ++i)
        r[i] = a[i] - b[i];
    return r;
}

inline double* prod(double* a, double* b, double* c) {
	double *r = new double[SIZE];
    for (uint i=0; i<SIZE; ++i)
        r[i] = a[i] * b[i] * c[i];
    return r;
}


inline double* vec_determinant(Armor::Array<double> &M) {
	return diff(
               sum( prod( M(0,0), M(1,1), M(2,2) ), prod( M(0,1), M(1,2), M(2,0) ), prod( M(0,2), M(1,0), M(2,1) ) ),
               sum( prod( M(2,0), M(1,1), M(0,2) ), prod( M(2,1), M(1,2), M(0,0) ), prod( M(2,2), M(1,0), M(0,1) ) )
			   );
}

#endif /* FEM_TOOLS_TEST_HH_ */
