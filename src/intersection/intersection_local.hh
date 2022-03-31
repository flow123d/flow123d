/*!
 *
﻿* Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    intersection_local.hh
 * @brief   Classes representing final intersection objects.
 * @author  Viktor Fris, Pavel Exner
 *
 */

#ifndef INTERSECTION_LOCAL_H_
#define INTERSECTION_LOCAL_H_

#include <armadillo>
#include "system/system.hh"
template <int spacedim> class ElementAccessor;

    
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
    /// Index of intersecting element in the component.
    unsigned int component_element_idx_;
    /// Index of intersecting element in the bulk.
    unsigned int bulk_element_idx_;

public:
    IntersectionLocalBase();
    /// Constructor taking in element indices.
    IntersectionLocalBase(unsigned int component_element_idx,
                          unsigned int bulk_element_idx);
    ~IntersectionLocalBase();
    
    unsigned int component_ele_idx() const; ///< Returns index of component element.
    unsigned int bulk_ele_idx() const;      ///< Returns index of bulk element.

    virtual double compute_measure() const =0;
};

/// First = element index, Second = pointer to intersection object.
typedef std::pair<unsigned int, IntersectionLocalBase*> ILpair;



inline unsigned int IntersectionLocalBase::component_ele_idx() const
{   return component_element_idx_; }

inline unsigned int IntersectionLocalBase::bulk_ele_idx() const
{   return bulk_element_idx_; }



/** @brief Class represents intersection of two elements.
 * 
 * It contains indices of intersecting elements (inherited from base class)
 * and vector of intersection points which provides barycentric coordinates
 * on both elements.
 */
template<unsigned int dimA, unsigned int dimB>
class IntersectionLocal : public IntersectionLocalBase
{
    /// Vector of intersection points.
    std::vector<IntersectionPoint<dimA,dimB>> i_points_;
    
public:

    /// Default constructor.
    IntersectionLocal();
    /// Constructor taking in element indices.
    IntersectionLocal(unsigned int component_element_idx, unsigned int bulk_element_idx);
    /// Copy constructor.
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
    double compute_measure() const override;
    

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
{   ASSERT(index < i_points_.size());
    return i_points_[index]; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionLocal<dimA,dimB>::size() const
{   return i_points_.size(); }




/** @brief Class represents an intersection point of simplex<N> and simplex<M>.
 * It contains barycentric coordinates of the point on both simplices.
 */
template<unsigned int dimA, unsigned int dimB>
class IntersectionPoint {
    
    arma::vec::fixed<dimA> comp_coords_; ///< Local coordinates of an IP on simplex<dimA>.
    arma::vec::fixed<dimB> bulk_coords_; ///< Local coordinates of an IP on simplex<dimB>.
    
public:

    IntersectionPoint();  ///< Default constructor.
    ~IntersectionPoint(); ///< Destructor.
    
    /// Constructs IP from the auxiliary one coming from computation.
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
    
    /**
     * Computes the real coordinates.
     * comp_ele is component element
     */

    arma::vec3 coords(ElementAccessor<3> comp_ele) const;
    
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



#endif /* INTERSECTION_LOCAL_H_ */
