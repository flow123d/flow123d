/*
 *  Created on: 14.3.2016
 *  Author: pex
 */

#include "intersectionaux.h"
#include "intersectionpoint.h"


namespace computeintersection{

template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::IntersectionAux(unsigned int component_element_idx,
                                            unsigned int bulk_element_idx,
                                            unsigned int comp_idx)
: component_element_idx_(component_element_idx), 
  bulk_element_idx_(bulk_element_idx),
  component_idx_(comp_idx)
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::IntersectionAux()
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::~IntersectionAux()
{}



// 1D-3D

template<>
double IntersectionAux<1,3>::compute_measure()
{
    //ASSERT(i_points_.size() > 1, "Not enough intersetion points to define a line.");
    double length = 0;
    
    if(i_points_.size() > 1)
    for(unsigned int i=0; i < i_points_.size()-1; i++)
    {
        length += std::abs(i_points_[i].local_bcoords_A()[1] - i_points_[i+1].local_bcoords_A()[1]);
    }
    return length;
}



// 2D-3D

template<>
double IntersectionAux<2,3>::compute_measure()
{
    double subtotal = 0.0;
    
    if(i_points_.size() > 2)
    for(unsigned int j = 2; j < i_points_.size();j++){
        //xprintf(Msg, "volani %d %d\n",j, i_points_.size());
        subtotal += fabs(i_points_[0].local_bcoords_A()(1)*(i_points_[j-1].local_bcoords_A()(2) - i_points_[j].local_bcoords_A()(2)) +
                 i_points_[j-1].local_bcoords_A()(1)*(i_points_[j].local_bcoords_A()(2) - i_points_[0].local_bcoords_A()(2)) +
                 i_points_[j].local_bcoords_A()(1)*(i_points_[0].local_bcoords_A()(2) - i_points_[j-1].local_bcoords_A()(2)));
    }
    return fabs(subtotal/2);
}


template<unsigned int dimA, unsigned int dimB>
ostream& operator<<(ostream& os, const IntersectionAux<dimA,dimB>& intersection)
{
    for(unsigned int i = 0; i < intersection.points().size(); i++)
        os << intersection.points()[i];
    
    return os;
}


template class IntersectionAux<1,2>;
template class IntersectionAux<2,2>;
template class IntersectionAux<1,3>;
template class IntersectionAux<2,3>;

template ostream& operator<< <1,2>(ostream &os, const IntersectionAux<1,2>& s);
template ostream& operator<< <2,2>(ostream &os, const IntersectionAux<2,2>& s);
template ostream& operator<< <1,3>(ostream &os, const IntersectionAux<1,3>& s);
template ostream& operator<< <2,3>(ostream &os, const IntersectionAux<2,3>& s);

}