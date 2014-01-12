/*
 * point.hh
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#ifndef POINT_HH_
#define POINT_HH_

#include <armadillo>


//#if HAVE_CXX11

/*
 * TODO:
 * need better resolution of various C++11 functionalities
 * e.g. following is supported from GCC 4.7
 */

//template <int spacedim>
//using Point = arma::vec::fixed<spacedim>

//#else

template<int spacedim>
class Space {
public:
    typedef typename arma::vec::fixed<spacedim> Point;
};
//#define

//#endif


#endif /* POINT_HH_ */
