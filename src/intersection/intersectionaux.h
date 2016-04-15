/*
 *  Created on: 14.3.2016
 *  Author: pex
 */

#ifndef INTERSECTIONAUX_H_
#define INTERSECTIONAUX_H_

#include "system/system.hh"

namespace computeintersection{

//forward declare
template<unsigned int, unsigned int> class IntersectionPointAux;
template<typename U, typename V> class ComputeIntersection;

/** @brief Internal auxiliary class representing intersecion object of of simplex<dimA> and simplex<dimB>.
 * 
 * It contains topology information and auxiliary intersection points.
 */
template<unsigned int dimA, unsigned int dimB>
class IntersectionAux{

    std::vector<IntersectionPointAux<dimA,dimB>> i_points_;
    unsigned int component_element_idx_;
    unsigned int bulk_element_idx_;
    bool pathologic_;
    
public:

    IntersectionAux();                                       ///< Default constructor.
    IntersectionAux(unsigned int component_element_idx, unsigned int bulk_element_idx);  ///< Constructor taking in element indices.
    virtual ~IntersectionAux();                                      ///< Destructor.

    /// Returns intersection points by a reference.
    std::vector<IntersectionPointAux<dimA,dimB>> &points();

    /// Returns intersection points by a constant reference.
    const std::vector<IntersectionPointAux<dimA,dimB>> &points() const;
    
    /// Returns intersection point of given @p index.
    const IntersectionPointAux<dimA,dimB> &operator[](unsigned int index) const;
    
    unsigned int size() const;              ///< Returns number of intersection points.
    unsigned int component_ele_idx() const; ///< Returns index of component element.
    unsigned int bulk_ele_idx() const;      ///< Returns index of bulk element.
    unsigned int is_pathologic() const;      ///< Returns index of bulk element.
    
    /// Computes the relative measure of intersection object.
    double compute_measure();
    
    /// Friend output operator.
    template<unsigned int dimAA, unsigned int dimBB>
    friend std::ostream& operator<<(std::ostream& os, const IntersectionAux<dimAA,dimBB>& intersection);
    
    template<typename U, typename V>
    friend class ComputeIntersection;
};

/********************************************* IMPLEMENTATION ***********************************************/

template<unsigned int dimA, unsigned int dimB>
inline std::vector< IntersectionPointAux< dimA, dimB > >& IntersectionAux<dimA,dimB>::points()
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const std::vector< IntersectionPointAux< dimA, dimB > >& IntersectionAux<dimA,dimB>::points() const
{   return i_points_; }

template<unsigned int dimA, unsigned int dimB>
inline const IntersectionPointAux< dimA, dimB >& IntersectionAux<dimA,dimB>::operator[](unsigned int index) const
{   ASSERT(index < i_points_.size(), "Index out of bounds.");
    return i_points_[index]; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionAux<dimA,dimB>::size() const
{   return i_points_.size(); }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionAux<dimA,dimB>::component_ele_idx() const
{   return component_element_idx_; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionAux<dimA,dimB>::bulk_ele_idx() const
{   return bulk_element_idx_; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionAux<dimA,dimB>::is_pathologic() const
{   return pathologic_; }


} // END NAMESPACE
#endif /* INTERSECTIONAUX_H_ */
