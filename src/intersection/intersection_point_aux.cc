/*
 * intersectionpoint.cpp
 *
 *  Created on: 11.4.2014
 *      Author: viktor
 */

// #include "mesh/elements.h"
#include "mesh/mesh.h"
#include "intersection_point_aux.hh"
#include "mesh/ref_element.hh"
// #include "mesh/elements.h"   //TODO what is the best way of include to use ElementFullIter ?
#include "mesh/mesh.h"
#include "system/system.hh"

using namespace std;


template<unsigned int N, unsigned int M>
void IntersectionPointAux<N,M>::clear()
{   local_bcoords_A_.zeros();
    local_bcoords_B_.zeros();
    idx_A_ = 0;
    idx_B_ = 0;
    orientation_ = IntersectionResult::none;
    dim_A_ = N;
    dim_B_ = M;
}

template<unsigned int N, unsigned int M>
IntersectionPointAux<N,M>::IntersectionPointAux()
{
    clear();
};

template<unsigned int N, unsigned int M>    
IntersectionPointAux<N,M>::IntersectionPointAux(const arma::vec::fixed<N+1> &lcA,
                                          const arma::vec::fixed<M+1> &lcB,
                                          unsigned int dim_A, 
                                          unsigned int dim_B)
    : local_bcoords_A_(lcA), local_bcoords_B_(lcB), dim_A_(dim_A), dim_B_(dim_B)
    {};


// template<unsigned int N, unsigned int M>
// IntersectionPointAux<N,M>::IntersectionPointAux(IntersectionPointAux<M, N> &IP){
//         local_bcoords_A_ = IP.local_bcoords_B();
//         local_bcoords_B_ = IP.local_bcoords_A();
//         idx_A_ = IP.idx_B();
//         idx_B_ = IP.idx_A();
//         orientation_ = IP.orientation();
//         dim_A_ = IP.dim_B();
//         dim_B_ = IP.dim_A();
//     };


template<unsigned int N, unsigned int M>
IntersectionPointAux<N,M>::IntersectionPointAux(const IntersectionPointAux<N,M-1> &IP, unsigned int idx_B){
    ASSERT_DBG(M>1 && M<4);
    
    local_bcoords_A_ = IP.local_bcoords_A();
    
    // FIXME: switch barycentric coords everywhere in intersections to remove this
    arma::vec::fixed<M> b;
    b.rows(0,M-2) = IP.local_bcoords_B().rows(1,M-1);
    b[M-1] = IP.local_bcoords_B()[0];
    
    local_bcoords_B_ = RefElement<M>::template interpolate<M-1>(b, idx_B);
    
    // switch back
    double bb = local_bcoords_B_[M];
    local_bcoords_B_.rows(1,M) = local_bcoords_B_.rows(0,M-1);
    local_bcoords_B_[0] = bb;
    
    dim_A_ = IP.dim_A();
    idx_A_ = IP.idx_A();
    orientation_ = IP.orientation();

    /**
     * TODO: set correct topology on B. Currently this is done ad hoc after call of this constructor.
     * Problem, dim_B_ can not be used as template parameter. Can we have some variant of interact without
     * template?
     */
    //dim_B_ = IP.dim_B();
    //idx_B_ =RefElement<M>::interact(Interaction<dim_B_, M-1>(IP.idx_B()));

    dim_B_ = M-1;
    idx_B_ = idx_B;
};


template<unsigned int N, unsigned int M>
IntersectionPointAux<N,M>::IntersectionPointAux(const IntersectionPointAux<N,M-2> &IP, unsigned int idx_B){
    ASSERT_DBG(M == 3);

    local_bcoords_A_ = IP.local_bcoords_A();
    
    // FIXME: switch barycentric coords everywhere in intersections to remove this
    arma::vec::fixed<2> b = IP.local_bcoords_B();
    std::swap(b[0],b[1]);
    
    local_bcoords_B_ = RefElement<3>::interpolate<1>(b, idx_B);
    
    // switch back
    double bb = local_bcoords_B_[M];
    local_bcoords_B_.rows(1,M) = local_bcoords_B_.rows(0,M-1);
    local_bcoords_B_[0] = bb;
    
    idx_A_ = IP.idx_A();
    idx_B_ = idx_B;
    orientation_ = IP.orientation();
    dim_A_ = IP.dim_A();
    dim_B_ = M-2;
};

template<unsigned int N, unsigned int M>
IntersectionPointAux<M,N> IntersectionPointAux<N,M>::switch_objects() const
{
    IntersectionPointAux<M,N> IP;
    IP.set_coordinates(local_bcoords_B_,local_bcoords_A_);
    IP.set_topology(idx_B_,dim_B_,idx_A_, dim_A_);
    IP.set_orientation(orientation_);
    return IP;
};


template<unsigned int N, unsigned int M>
arma::vec::fixed< 3  > IntersectionPointAux<N,M>::coords(ElementFullIter ele) const
{
    ASSERT_DBG(N == ele->dim());
    
    arma::vec::fixed< 3  > c;
    c.zeros();
    for(unsigned int i=0; i<N+1; i++)
        c += local_bcoords_A_[i]*ele->node[i]->point();
        
    return c;
}
 
/*
template<> bool IntersectionPointAux<2,3>::operator<(const IntersectionPointAux<2,3> &ip) const{
	return local_bcoords_A_[1] < ip.local_bcoords_A()[1] ||     // compare by x coordinate
           (local_bcoords_A_[1] == ip.local_bcoords_A()[1] &&   // in case of tie
           local_bcoords_A_[2] < ip.local_bcoords_A()[2]);      // compare by y coordinate
};
*/



template<unsigned int N, unsigned int M> ostream& operator<<(ostream& os, const IntersectionPointAux< N,M >& s)
{
    os << "Local coords on element A(id=" << s.idx_A_ << ", dim=" << s.dim_A_ << ")" << endl;
    s.local_bcoords_A_.print(os);
    os << "Local coords on element B(id=" << s.idx_B_ << ", dim=" << s.dim_B_ << ")" << endl;
    s.local_bcoords_B_.print(os);
    os << "Orientation: " << int(s.orientation_) << " Patological: " << s.is_pathologic() << endl;
    return os;
}

template class IntersectionPointAux<1,2>;
template class IntersectionPointAux<2,1>;
template class IntersectionPointAux<2,2>;
template class IntersectionPointAux<1,3>;
template class IntersectionPointAux<3,1>;
template class IntersectionPointAux<2,3>;
template class IntersectionPointAux<3,2>;

template ostream& operator<< <1,2>(ostream &os, const IntersectionPointAux<1,2>& s);
template ostream& operator<< <2,1>(ostream &os, const IntersectionPointAux<2,1>& s);
template ostream& operator<< <2,2>(ostream &os, const IntersectionPointAux<2,2>& s);
template ostream& operator<< <1,3>(ostream &os, const IntersectionPointAux<1,3>& s);
template ostream& operator<< <3,1>(ostream &os, const IntersectionPointAux<3,1>& s);
template ostream& operator<< <2,3>(ostream &os, const IntersectionPointAux<2,3>& s);
template ostream& operator<< <3,2>(ostream &os, const IntersectionPointAux<3,2>& s);





