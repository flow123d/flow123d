/*
 * intersectionpoint.cpp
 *
 *  Created on: 11.4.2014
 *      Author: viktor
 */

#include "intersectionpoint.h"

namespace computeintersection{


/*template<> arma::vec::fixed<4> interpolate<1,3>(arma::vec::fixed<2> &coord, unsigned int sub_simplex_idx){

   // std::array<vec::fixed<4>, 2 >  simplex_M_vertices = RefElement<3>.sub_element<1>.bary_coords(sub_simplex_idx);
    arma::vec::fixed<4> sum;
    sum.zeros();
   // for(int i=0; i<2; i++) sum += coord[i]*simplex_M_vertices(i);
    return sum;
};*/


/*template<int M, int N> arma::vec::fixed<N+1> interpolate(arma::vec::fixed<M+1> &coord, unsigned int sub_simplex_idx){

	std::array<arma::vec::fixed<N+1>, M+1> simplex_M_vertices = RefSimplex<N>::bary_coords<M>(sub_simplex_idx);
	    arma::vec::fixed<N+1> sum;
	    sum.zeros();
	    for(int i=0; i<M+1; i++) sum += coord[i]*simplex_M_vertices[i];
	    return sum;
};*/

} // END namespace



