/*
 * intersectionpoint.cpp
 *
 *  Created on: 11.4.2014
 *      Author: viktor
 */

#include "intersectionpoint.h"
#include "mesh/ref_element.hh"
#include "system/system.hh"

using namespace std;
namespace computeintersection{

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::clear()
{   local_coords1_.zeros();
    local_coords2_.zeros();
    side_idx1_ = 0;
    side_idx2_ = 0;
    orientation_ = 1;
    pathologic_ = false;
    dim_A_ = N;
    dim_B_ = M;
}

template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint()
{
    clear();
};

template<unsigned int N, unsigned int M>    
IntersectionPoint<N,M>::IntersectionPoint(const arma::vec::fixed<N+1> &lc1,
                                          const arma::vec::fixed<M+1> &lc2,
                                          unsigned int dim_A, 
                                          unsigned int dim_B)
    : local_coords1_(lc1), local_coords2_(lc2), dim_A_(dim_A), dim_B_(dim_B)
    {};


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<M, N> &IP){
        local_coords1_ = IP.local_coords2();
        local_coords2_ = IP.local_coords1();
        side_idx1_ = IP.side_idx2();
        side_idx2_ = IP.side_idx1();
        orientation_ = IP.orientation();
        pathologic_ = IP.is_pathologic();
        dim_A_ = IP.dim_B();
        dim_B_ = IP.dim_A();
    };


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-1> &IP, unsigned int side_idx2){
    ASSERT(M>1 && M<4,"Wrong the second dimension in an IntersectionPoint (allowed 2 and 3 only)");
    
    local_coords1_ = IP.local_coords1();
    local_coords2_ = RefElement<M>::template interpolate<M-1>(IP.local_coords2(), side_idx2);
    side_idx1_ = IP.side_idx1();
    side_idx2_ = side_idx2;
    orientation_ = IP.orientation();
    pathologic_ = IP.is_pathologic();
    dim_A_ = IP.dim_A();
    dim_B_ = M-1;
};


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-2> &IP, unsigned int side_idx2){
    ASSERT(M == 3,"Wrong the second dimension in an IntersectionPoint (allowed 3 only)");

    local_coords1_ = IP.local_coords1();
    local_coords2_ = RefElement<3>::interpolate<1>(IP.local_coords2(), side_idx2);
    side_idx1_ = IP.side_idx1();
    side_idx2_ = side_idx2;
    orientation_ = IP.orientation();
    pathologic_ = IP.is_pathologic();
    dim_A_ = IP.dim_A();
    dim_B_ = M-2;
};

    
template<> bool IntersectionPoint<2,3>::operator<(const IntersectionPoint<2,3> &ip) const{
	return local_coords1_[1] < ip.local_coords1()[1] ||
		(fabs((double)(local_coords1_[1] - ip.local_coords1()[1])) < 0.00000001 &&
		 local_coords1_[2] < ip.local_coords1()[2]);
};

template<unsigned int N, unsigned int M> ostream& operator<<(ostream& os, const IntersectionPoint< N,M >& s)
{
    os << "Local coords on element A(id=" << s.side_idx1_ << ", dim=" << s.dim_A_ << ")" << endl;
    s.local_coords1_.print(os);
    os << "Local coords on element B(id=" << s.side_idx2_ << ", dim=" << s.dim_B_ << ")" << endl;
    s.local_coords2_.print(os);
    os << "Orientation: " << s.orientation_ << " Patological: " << s.pathologic_ << endl;
    return os;
}

template class IntersectionPoint<1,2>;
template class IntersectionPoint<1,3>;
template class IntersectionPoint<2,1>;
template class IntersectionPoint<2,3>;
template class IntersectionPoint<3,1>;
template class IntersectionPoint<3,2>;

template ostream& operator<< <1,2>(ostream &os, const IntersectionPoint<1,2>& s); 
template ostream& operator<< <1,3>(ostream &os, const IntersectionPoint<1,3>& s); 
template ostream& operator<< <2,1>(ostream &os, const IntersectionPoint<2,1>& s); 
template ostream& operator<< <2,3>(ostream &os, const IntersectionPoint<2,3>& s); 
template ostream& operator<< <3,1>(ostream &os, const IntersectionPoint<3,1>& s); 
template ostream& operator<< <3,2>(ostream &os, const IntersectionPoint<3,2>& s); 

} // END namespace



