/*
 *  Created on: 14.3.2016
 *  Author: pex
 */

#include "intersection_aux.hh"
#include "intersection_point_aux.hh"
#include "mesh/ref_element.hh"


template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::IntersectionAux(unsigned int component_element_idx,
                                            unsigned int bulk_element_idx)
: component_element_idx_(component_element_idx), 
  bulk_element_idx_(bulk_element_idx),
  n_duplicities_(0)
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::IntersectionAux()
{}

template<unsigned int dimA, unsigned int dimB>
IntersectionAux<dimA,dimB>::~IntersectionAux()
{}


template<unsigned int dimA, unsigned int dimB>
unsigned int IntersectionAux<dimA,dimB>::ips_on_single_object() const
{
    const uint invalid_face = -1;
    if(size() < dimB) return invalid_face;
    
    ASSERT_DBG((dimA == 1 && dimB == 2) || (dimA == 2 && dimB == 3));
    
    //test if all IPs lie in a single face of tetrahedron
    vector<unsigned int> face_counter(RefElement<dimB>::n_sides, 0);
    
    for(const IntersectionPointAux<dimA,dimB>& ipf : i_points_) {
        switch(ipf.dim_B()){
            case 0: {
                for(uint i=0; i<RefElement<dimB>::n_sides_per_node; i++)
                    face_counter[RefElement<dimB>::interact(Interaction<dimA,0>(ipf.idx_B()))[i]]++;
                break;
            }
            case 1: {
                if(dimB == 3)
                    for(uint i=0; i<RefElement<dimB>::n_sides_per_line; i++)
                        face_counter[RefElement<dimB>::interact(Interaction<dimA,1>(ipf.idx_B()))[i]]++;
                else
                    face_counter[ipf.idx_B()]++;
                break;
            }
            case 2:
                if(dimB == 3)
                    face_counter[ipf.idx_B()]++;
                else 
                    return invalid_face;    // cannot be on face
                break;
            case 3: 
                return invalid_face;    // cannot be on face
                break;
            default: ASSERT(0)(ipf.dim_B());
        }
    }
    
    for(uint f=0; f< RefElement<dimB>::n_sides; f++){
        if(face_counter[f] == size()){ // all IPs lie in a single face
            DBGCOUT(<< "all IPs in face: " << f <<"\n");
            return f;
        }
    }
    //should never reach
    return invalid_face;
}


// 1D-3D

template<>
double IntersectionAux<1,3>::compute_measure()
{
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


