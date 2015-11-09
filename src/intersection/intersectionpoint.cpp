/*
 * intersectionpoint.cpp
 *
 *  Created on: 11.4.2014
 *      Author: viktor
 */

#include "intersectionpoint.h"
#include "mesh/ref_element.hh"

namespace computeintersection{

template<int N, int M>
IntersectionPoint<N,M>::IntersectionPoint()
{
    clear();
};

template<int N, int M>    
IntersectionPoint<N,M>::IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
        const arma::vec::fixed<M+1> &lc2,
        int side1,
        int side2,
        unsigned int ori,
        bool vertex,
        bool patological)
        : local_coords1(lc1),
          local_coords2(lc2),
          side_idx1(side1),
          side_idx2(side2),
          orientation(ori),
          is_vertex_(vertex),
          is_patological_(patological) {};


template<int N, int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<M, N> &IP){
        local_coords1 = IP.get_local_coords2();
        local_coords2 = IP.get_local_coords1();
        side_idx1 = IP.get_side2();
        side_idx2 = IP.get_side1();
        orientation = IP.get_orientation();
        is_vertex_ = IP.is_vertex();
        is_patological_ = IP.is_patological();
    };


template<int N, int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-1> &IP){
    ASSERT(M>1 && M<4,"Wrong the second dimension in an IntersectionPoint (allowed 2 and 3 only)");
    
    local_coords1 = IP.get_local_coords1();
    local_coords2 = RefElement<M>::template interpolate<M-1>(IP.get_local_coords2(), IP.get_side2());
    side_idx1 = IP.get_side1();
    side_idx2 = IP.get_side2();
    orientation = IP.get_orientation();
    is_vertex_ = IP.is_vertex();
    is_patological_ = IP.is_patological();
};


template<int N, int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-2> &IP){
    ASSERT(M == 3,"Wrong the second dimension in an IntersectionPoint (allowed 3 only)");

    local_coords1 = IP.get_local_coords1();
    local_coords2 = RefElement<3>::interpolate<1>(IP.get_local_coords2(), IP.get_side2());
    side_idx1 = IP.get_side1();
    side_idx2 = IP.get_side2();
    orientation = IP.get_orientation();
    is_vertex_ = IP.is_vertex();
    is_patological_ = IP.is_patological();
};

    
template<> bool IntersectionPoint<2,3>::operator<(const IntersectionPoint<2,3> &ip) const{
	return local_coords1[1] < ip.get_local_coords1()[1] ||
		(fabs((double)(local_coords1[1] - ip.get_local_coords1()[1])) < 0.00000001 &&
		 local_coords1[2] < ip.get_local_coords1()[2]);
};


template class IntersectionPoint<1,2>;
template class IntersectionPoint<1,3>;
template class IntersectionPoint<2,1>;
template class IntersectionPoint<2,3>;
template class IntersectionPoint<3,1>;
template class IntersectionPoint<3,2>;

} // END namespace



