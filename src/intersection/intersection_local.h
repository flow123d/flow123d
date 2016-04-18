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

    
//forwward declare
template<unsigned int, unsigned int> class IntersectionPointAux;
template<unsigned int, unsigned int> class IntersectionAux;

template<unsigned int, unsigned int> class IntersectionPoint;
template<unsigned int, unsigned int> class IntersectionLocal;
template<unsigned int dimA, unsigned int dimB> std::ostream& operator<<(std::ostream& os, const IntersectionLocal<dimA,dimB>& il);
template<unsigned int dimA, unsigned int dimB> std::ostream& operator<<(std::ostream& os, const IntersectionPoint<dimA,dimB>& ip);


/** @brief Common base for intersection object.
 * 
 * This base class provides unification of all intersection objects.
 * The common part is a pair of component element index and bulk element index.
 * 
 * The derived class @p IntersectionLocal<dimA,dimB> differs in coordinates length
 * according to @p dimA and @p dimB.
 * 
 */
class IntersectionLocalBase
{
protected:
    unsigned int component_element_idx_;
    unsigned int bulk_element_idx_;
    unsigned int component_idx_;

public:
    IntersectionLocalBase(){}
    /// Constructor taking in element indices.
    IntersectionLocalBase(unsigned int component_element_idx,
                          unsigned int bulk_element_idx,
                          unsigned int component_idx);
    ~IntersectionLocalBase(){}
    
    unsigned int component_ele_idx() const; ///< Returns index of component element.
    unsigned int bulk_ele_idx() const;      ///< Returns index of bulk element.
    unsigned int component_idx() const;     ///< Returns index of component.
};


inline unsigned int IntersectionLocalBase::component_ele_idx() const
{   return component_element_idx_; }

inline unsigned int IntersectionLocalBase::bulk_ele_idx() const
{   return bulk_element_idx_; }

inline unsigned int IntersectionLocalBase::component_idx() const
{   return component_idx_; }



/** @brief Class represents intersection of two elements.
 * 
 * It contains indices of intersecting elements (inherited from base class)
 * and vector of intersection points which provides barycentric coordinates
 * on both elements.
 */
template<unsigned int dimA, unsigned int dimB>
class IntersectionLocal : public IntersectionLocalBase
{
    /// Vector of intersectio points.
    std::vector<IntersectionPoint<dimA,dimB>> i_points_;
    
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
    std::vector<IntersectionPoint<dimA,dimB>> &points();

    /// Returns intersection points by a constant reference.
    const std::vector<IntersectionPoint<dimA,dimB>> &points() const;

    /// Returns intersection point of given @p index.
    const IntersectionPoint<dimA,dimB> &operator[](unsigned int index) const;
    
    unsigned int size() const;              ///< Returns number of intersection points.
    //@}
    
    /// Computes the relative measure of intersection object.
    double compute_measure();
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionLocal<dimA,dimB>& intersection);
};

/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int dimA, unsigned int dimB>
inline std::vector< IntersectionPoint< dimA, dimB > >& IntersectionLocal<dimA,dimB>::points()
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const std::vector< IntersectionPoint< dimA, dimB > >& IntersectionLocal<dimA,dimB>::points() const
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const IntersectionPoint< dimA, dimB >& IntersectionLocal<dimA,dimB>::operator[](unsigned int index) const
{   ASSERT(index < i_points_.size(), "Index out of bounds.");
    return i_points_[index]; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionLocal<dimA,dimB>::size() const
{   return i_points_.size(); }




/**
 * Class represents an intersection point of simplex<N> and simplex<M>.
 * It contains barycentric coordinates of the point on both simplices.
 */
template<unsigned int dimA, unsigned int dimB> class IntersectionPoint {
    
    arma::vec::fixed<dimA> comp_coords_; ///< Local coordinates of an IP on simplex<dimA>.
    arma::vec::fixed<dimB> bulk_coords_; ///< Local coordinates of an IP on simplex<dimB>.
    
public:

    IntersectionPoint();  ///< Default constructor.
    ~IntersectionPoint(); ///< Destructor.
    
    IntersectionPoint(const IntersectionPointAux<dimA,dimB> &p);
    /**
     * Constructor taking local coordinates on simplices as input parameters.
     * @param comp_coords local coordinates of IP in Simplex<dimA>
     * @param bulk_coords local coordinates of IP in Simplex<dimB>
     */
    IntersectionPoint(const arma::vec::fixed<dimA> &comp_coords, const arma::vec::fixed<dimB> &bulk_coords);
    
    ///@name Getters.
    //@{
    /// Returns local coordinates in the Simplex<N>.
    const arma::vec::fixed<dimA> &comp_coords() const;
    
    /// Returns local coordinates in the Simplex<M>.
    const arma::vec::fixed<dimB> &bulk_coords() const;
    //@}
    
    /// Computes the real coordinates.
    arma::vec3 coords(ElementFullIter comp_ele) const;
    
    /// Friend output operator.
    friend std::ostream& operator<< <>(std::ostream& os, const IntersectionPoint<dimA,dimB>& IP);
};


/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int dimA, unsigned int dimB>
const arma::vec::fixed< dimA  >& IntersectionPoint<dimA,dimB>::comp_coords() const
{   return comp_coords_; }

template<unsigned int dimA, unsigned int dimB>
const arma::vec::fixed< dimB  >& IntersectionPoint<dimA,dimB>::bulk_coords() const
{   return bulk_coords_; }

} // END NAMESPACE
#endif /* INTERSECTION_LOCAL_H_ */
