/*
 * intersectionpoint.cpp
 *
 *  Created on: 11.4.2014
 *      Author: viktor
 */

#include "intersectionpoint.h"
#include "mesh/ref_element.hh"
#include <mesh/elements.h>
#include "system/system.hh"

using namespace std;
namespace computeintersection{

template<unsigned int N, unsigned int M>
void IntersectionPoint<N,M>::clear()
{   local_bcoords_A_.zeros();
    local_bcoords_B_.zeros();
    idx_A_ = 0;
    idx_B_ = 0;
    orientation_ = 1;
    dim_A_ = N;
    dim_B_ = M;
}

template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint()
{
    clear();
};

template<unsigned int N, unsigned int M>    
IntersectionPoint<N,M>::IntersectionPoint(const arma::vec::fixed<N+1> &lcA,
                                          const arma::vec::fixed<M+1> &lcB,
                                          unsigned int dim_A, 
                                          unsigned int dim_B)
    : local_bcoords_A_(lcA), local_bcoords_B_(lcB), dim_A_(dim_A), dim_B_(dim_B)
    {};


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<M, N> &IP){
        local_bcoords_A_ = IP.local_bcoords_B();
        local_bcoords_B_ = IP.local_bcoords_A();
        idx_A_ = IP.idx_B();
        idx_B_ = IP.idx_A();
        orientation_ = IP.orientation();
        dim_A_ = IP.dim_B();
        dim_B_ = IP.dim_A();
    };


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-1> &IP, unsigned int idx_B){
    ASSERT(M>1 && M<4,"Wrong the second dimension in an IntersectionPoint (allowed 2 and 3 only)");
    
    local_bcoords_A_ = IP.local_bcoords_A();
    local_bcoords_B_ = RefElement<M>::template interpolate<M-1>(IP.local_bcoords_B(), idx_B);
    idx_A_ = IP.idx_A();
    idx_B_ = idx_B;
    orientation_ = IP.orientation();
    dim_A_ = IP.dim_A();
    dim_B_ = M-1;
};


template<unsigned int N, unsigned int M>
IntersectionPoint<N,M>::IntersectionPoint(IntersectionPoint<N,M-2> &IP, unsigned int idx_B){
    ASSERT(M == 3,"Wrong the second dimension in an IntersectionPoint (allowed 3 only)");

    local_bcoords_A_ = IP.local_bcoords_A();
    local_bcoords_B_ = RefElement<3>::interpolate<1>(IP.local_bcoords_B(), idx_B);
    idx_A_ = IP.idx_A();
    idx_B_ = idx_B;
    orientation_ = IP.orientation();
    dim_A_ = IP.dim_A();
    dim_B_ = M-2;
};

template<unsigned int N, unsigned int M>
arma::vec::fixed< 3  > IntersectionPoint<N,M>::coords(ElementFullIter ele)
{
    ASSERT(N == ele->dim(), "Element vs intersection point dimension mismatch.");
    
    arma::vec::fixed< 3  > c;
    c.zeros();
    for(unsigned int i=0; i<N+1; i++)
        c += local_bcoords_A_[i+1]*ele->node[i]->point();
        
    return c;
}
 

template<> bool IntersectionPoint<2,3>::operator<(const IntersectionPoint<2,3> &ip) const{
	return local_bcoords_A_[1] < ip.local_bcoords_A()[1] ||     // compare by x coordinate
           (local_bcoords_A_[1] == ip.local_bcoords_A()[1] &&   // in case of tie
           local_bcoords_A_[2] < ip.local_bcoords_A()[2]);      // compare by y coordinate
};

template<unsigned int N, unsigned int M> ostream& operator<<(ostream& os, const IntersectionPoint< N,M >& s)
{
    os << "Local coords on element A(id=" << s.idx_A_ << ", dim=" << s.dim_A_ << ")" << endl;
    s.local_bcoords_A_.print(os);
    os << "Local coords on element B(id=" << s.idx_B_ << ", dim=" << s.dim_B_ << ")" << endl;
    s.local_bcoords_B_.print(os);
    os << "Orientation: " << s.orientation_ << " Patological: " << s.is_pathologic() << endl;
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



