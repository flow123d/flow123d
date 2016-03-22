/*
 *  Created on: 14.3.2016
 *  Author: pex
 */

#ifndef INTERSECTION_LOCAL_H_
#define INTERSECTION_LOCAL_H_

#include <armadillo>
#include <mesh/mesh_types.hh>
#include "system/system.hh"

namespace computeintersection{


class IntersectionLocalBase
{
protected:
    unsigned int component_element_idx_;
    unsigned int bulk_element_idx_;

public:
    IntersectionLocalBase(){}
    /// Constructor taking in element indices.
    IntersectionLocalBase(unsigned int component_element_idx, unsigned int bulk_element_idx);
    ~IntersectionLocalBase(){}
    
    unsigned int component_ele_idx() const; ///< Returns index of component element.
    unsigned int bulk_ele_idx() const;      ///< Returns index of bulk element.
};


inline unsigned int IntersectionLocalBase::component_ele_idx() const
{   return component_element_idx_; }

inline unsigned int IntersectionLocalBase::bulk_ele_idx() const
{   return bulk_element_idx_; }



//forwward declare
template<unsigned int, unsigned int> class IntersectionPoint;
template<unsigned int, unsigned int> class IntersectionAux;

template<unsigned int, unsigned int> class IntersectionPointX;
template<unsigned int, unsigned int> class IntersectionLocal;
template<unsigned int dimA, unsigned int dimB> std::ostream& operator<<(std::ostream& os, const IntersectionLocal<dimA,dimB>& il);
template<unsigned int dimA, unsigned int dimB> std::ostream& operator<<(std::ostream& os, const IntersectionPointX<dimA,dimB>& ip);

template<unsigned int dimA, unsigned int dimB>
class IntersectionLocal : public IntersectionLocalBase
{
    std::vector<IntersectionPointX<dimA,dimB>> i_points_;
    
public:

    /// Default constructor.
    IntersectionLocal();
    /// Constructor taking in element indices.
    IntersectionLocal(unsigned int component_element_idx, unsigned int bulk_element_idx);
    IntersectionLocal(const IntersectionAux<dimA, dimB> &iaux);
    /// Destructor.
    ~IntersectionLocal();

    ///@name Getters.
    //@{
    /// Returns intersection points by a reference.
    std::vector<IntersectionPointX<dimA,dimB>> &points();

    /// Returns intersection points by a constant reference.
    const std::vector<IntersectionPointX<dimA,dimB>> &points() const;

    /// Returns intersection point of given @p index.
    const IntersectionPointX<dimA,dimB> &operator[](unsigned int index) const;
    
    unsigned int size() const;              ///< Returns number of intersection points.
    //@}
    
    /// Computes the relative measure of intersection object.
    double compute_measure();
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionLocal<dimA,dimB>& intersection);
};

/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int dimA, unsigned int dimB>
inline std::vector< IntersectionPointX< dimA, dimB > >& IntersectionLocal<dimA,dimB>::points()
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const std::vector< IntersectionPointX< dimA, dimB > >& IntersectionLocal<dimA,dimB>::points() const
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const IntersectionPointX< dimA, dimB >& IntersectionLocal<dimA,dimB>::operator[](unsigned int index) const
{   ASSERT(index < i_points_.size(), "Index out of bounds.");
    return i_points_[index]; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionLocal<dimA,dimB>::size() const
{   return i_points_.size(); }




/**
 * Class represents an intersection point of simplex<N> and simplex<M>.
 * It contains barycentric coordinates of the point on both simplices.
 */
template<unsigned int dimA, unsigned int dimB> class IntersectionPointX {
    
    arma::vec::fixed<dimA+1> comp_bcoords_; ///< Barycentric coordinates of an IP on simplex<N>.
    arma::vec::fixed<dimB+1> bulk_bcoords_; ///< Barycentric coordinates of an IP on simplex<M>.
    
public:

    IntersectionPointX(){};  ///< Default constructor.
    ~IntersectionPointX(){}; ///< Destructor.
    
    IntersectionPointX(const IntersectionPoint<dimA,dimB> &p);
    /**
     * Constructor taking barycentric coordinates on simplices as input parameters.
     * @param comp_bcoords barycentric coordinates of IP in Simplex<dimA>
     * @param bulk_bcoords barycentric coordinates of IP in Simplex<dimB>
     */
    IntersectionPointX(const arma::vec::fixed<dimA+1> &comp_bcoords, const arma::vec::fixed<dimB+1> &bulk_bcoords);
    
    ///@name Getters.
    //@{
    /// Returns barycentric coordinates in the Simplex<N>.
    const arma::vec::fixed<dimA+1> &comp_bcoords() const;
    
    /// Returns barycentric coordinates in the Simplex<M>.
    const arma::vec::fixed<dimB+1> &bulk_bcoords() const;
    //@}
    
    
    arma::vec3 coords(ElementFullIter comp_ele) const;
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionPointX<dimA,dimB>& IP);
};


/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int dimA, unsigned int dimB>
const arma::vec::fixed< dimA + 1  >& IntersectionPointX<dimA,dimB>::comp_bcoords() const
{   return comp_bcoords_; }

template<unsigned int dimA, unsigned int dimB>
const arma::vec::fixed< dimB + 1  >& IntersectionPointX<dimA,dimB>::bulk_bcoords() const
{   return bulk_bcoords_; }

} // END NAMESPACE
#endif /* INTERSECTION_LOCAL_H_ */
