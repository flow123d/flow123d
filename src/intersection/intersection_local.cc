
#include "intersection_local.hh"
#include "intersection_aux.hh"
#include "intersection_point_aux.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"

#include <iostream>

using namespace std;

IntersectionLocalBase::IntersectionLocalBase()
{}


IntersectionLocalBase::~IntersectionLocalBase()
{}

IntersectionLocalBase::IntersectionLocalBase(unsigned int component_element_idx, 
                                             unsigned int bulk_element_idx)
: component_element_idx_(component_element_idx), bulk_element_idx_(bulk_element_idx)
{}


template<unsigned int dimA, unsigned int dimB>
IntersectionLocal<dimA,dimB>::IntersectionLocal()
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionLocal<dimA,dimB>::IntersectionLocal(unsigned int component_element_idx, unsigned int bulk_element_idx)
: IntersectionLocalBase(component_element_idx, bulk_element_idx)
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionLocal<dimA,dimB>::IntersectionLocal(const IntersectionAux< dimA, dimB >& iaux)
: IntersectionLocalBase(iaux.component_ele_idx(), iaux.bulk_ele_idx())
{
    i_points_.resize(iaux.size());
    for(unsigned int i = 0; i < iaux.size(); i++)
    {
        i_points_[i] = IntersectionPoint<dimA,dimB>(iaux[i]);
    }
}


template<unsigned int dimA, unsigned int dimB>
IntersectionLocal<dimA,dimB>::~IntersectionLocal()
{}


template<unsigned int dimA, unsigned int dimB>
double IntersectionLocal<dimA,dimB>::compute_measure() const
{
    ASSERT_LT_DBG(i_points_.size(), 3 ); // avoid 2d-3d case and degenerated 2d-2d
    double length = 0;

    if(i_points_.size() > 1) // zero measure for point intersections
        for(unsigned int i=0; i < i_points_.size()-1; i++)
        {
            length += abs(i_points_[i].comp_coords()[0] - i_points_[i+1].comp_coords()[0]);
        }
    return length;
}



// 2D-3D
template<>
double IntersectionLocal<2,3>::compute_measure() const
{
    double subtotal = 0.0;
    
    if(i_points_.size() > 2)
    for(unsigned int j = 2; j < i_points_.size();j++){
        //DebugOut().fmt("volani {} {}\n", j, i_points_.size());
        subtotal += fabs(i_points_[0].comp_coords()(0)*(i_points_[j-1].comp_coords()(1) - i_points_[j].comp_coords()(1)) +
                 i_points_[j-1].comp_coords()(0)*(i_points_[j].comp_coords()(1) - i_points_[0].comp_coords()(1)) +
                 i_points_[j].comp_coords()(0)*(i_points_[0].comp_coords()(1) - i_points_[j-1].comp_coords()(1)));
    }
    return fabs(subtotal/2);
}



template<unsigned int dimA, unsigned int dimB>
IntersectionPoint<dimA,dimB>::IntersectionPoint()
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionPoint<dimA,dimB>::~IntersectionPoint()
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionPoint<dimA,dimB>::IntersectionPoint(const IntersectionPointAux<dimA,dimB>& p)
: comp_coords_(p.local_bcoords_A().subvec(1,dimA)), bulk_coords_(p.local_bcoords_B().subvec(1,dimB))
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionPoint<dimA,dimB>::IntersectionPoint(const arma::vec::fixed< dimA  >& comp_coords, 
                                                const arma::vec::fixed< dimB  >& bulk_coords)
: comp_coords_(comp_coords), bulk_coords_(bulk_coords)
{}



template<unsigned int dimA, unsigned int dimB>
arma::vec3 IntersectionPoint<dimA,dimB>::coords(ElementAccessor<3> comp_ele) const
{
    ASSERT_DBG(dimA == comp_ele.dim());
    
    arma::vec3 c;
    c.zeros();
    double complement = 1.0;
    for(unsigned int i=0; i<dimA; i++)
    {
        c += comp_coords_[i] * (*comp_ele.node(i+1));
        complement -= comp_coords_[i];
    }
    c += complement * (*comp_ele.node(0));
        
    return c;
}




template<unsigned int dimA, unsigned int dimB> ostream& operator<<(ostream& os, const IntersectionLocal<dimA,dimB>& il)
{
    os << "IntersectionLocal<" << dimA << "," << dimB << ">: c " << il.component_element_idx_ << ", b " << 
        il.bulk_element_idx_ << ", size " << il.i_points_.size() << endl;
    for (unsigned int i = 0; i < il.i_points_.size(); i++) {
        os << i << ": " << il[i] << endl;
    }        
    return os;
}

template<unsigned int dimA, unsigned int dimB> ostream& operator<<(ostream& os, const IntersectionPoint<dimA,dimB>& ip)
{
    os << "[";
    for(unsigned j= 0; j < dimA-1; j++)
        os << ip.comp_coords_[j] << " ";
    os << ip.comp_coords_[dimA-1] << "]\t[";
    
    for(unsigned j= 0; j < dimB-1; j++)
        os << ip.bulk_coords_[j] << " ";
    os << ip.bulk_coords_[dimB-1] << "]";
    os << endl;
    
    return os;
}


template class IntersectionPoint<1,2>;
template class IntersectionPoint<2,2>;
template class IntersectionPoint<1,3>;
template class IntersectionPoint<2,3>;

template class IntersectionLocal<1,2>;
template class IntersectionLocal<2,2>;
template class IntersectionLocal<1,3>;
template class IntersectionLocal<2,3>;

template ostream& operator<< <1,2>(ostream &os, const IntersectionPoint<1,2>& s);
template ostream& operator<< <2,2>(ostream &os, const IntersectionPoint<2,2>& s);
template ostream& operator<< <1,3>(ostream &os, const IntersectionPoint<1,3>& s); 
template ostream& operator<< <2,3>(ostream &os, const IntersectionPoint<2,3>& s);

template ostream& operator<< <1,2>(ostream &os, const IntersectionLocal<1,2>& s);
template ostream& operator<< <2,2>(ostream &os, const IntersectionLocal<2,2>& s);
template ostream& operator<< <1,3>(ostream &os, const IntersectionLocal<1,3>& s); 
template ostream& operator<< <2,3>(ostream &os, const IntersectionLocal<2,3>& s); 
