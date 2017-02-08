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
 * @file    intersection_aux.hh
 * @brief   Internal class representing intersection object.
 * @author  Pavel Exner
 *
 */

#ifndef INTERSECTIONAUX_H_
#define INTERSECTIONAUX_H_

#include "system/system.hh"

namespace computeintersection{

//forward declare
template<unsigned int, unsigned int> class IntersectionPointAux;
template<typename U, typename V> class ComputeIntersection;

/** @brief Internal auxiliary class representing intersection object of simplex<dimA> and simplex<dimB>.
 * 
 * It contains topology information and auxiliary intersection points.
 * Used in ComputeIntersection classes.
 */
template<unsigned int dimA, unsigned int dimB>
class IntersectionAux{

    /// Vector of internal intersection points.
    std::vector<IntersectionPointAux<dimA,dimB>> i_points_;

    /// Index of intersecting element in the component.
    unsigned int component_element_idx_;
    /// Index of intersecting element in the bulk.
    unsigned int bulk_element_idx_;
    /// Index of the intersecting component.
    unsigned int component_idx_;
    /// Flag for pathologic case.
    bool pathologic_;
    
public:

    /// Default constructor.
    IntersectionAux();
    /// Constructor taking in element indices.
    IntersectionAux(unsigned int component_element_idx,
                    unsigned int bulk_element_idx);
    /// Destructor.
    virtual ~IntersectionAux();

    /// Returns intersection points by a reference.
    std::vector<IntersectionPointAux<dimA,dimB>> &points();

    /// Returns intersection points by a constant reference.
    const std::vector<IntersectionPointAux<dimA,dimB>> &points() const;
    
    /// Returns intersection point of given @p index.
    const IntersectionPointAux<dimA,dimB> &operator[](unsigned int index) const;
    
    void set_component_idx(unsigned int comp_idx);
    
    unsigned int size() const;              ///< Returns number of intersection points.
    unsigned int component_ele_idx() const; ///< Returns index of component element.
    unsigned int bulk_ele_idx() const;      ///< Returns index of bulk element.
    unsigned int component_idx() const;     ///< Returns index of component.
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
{   ASSERT_DBG(index < i_points_.size());
    return i_points_[index]; }

template<unsigned int dimA, unsigned int dimB>
inline void IntersectionAux<dimA,dimB>::set_component_idx(unsigned int comp_idx)
{   component_idx_ = comp_idx; }

    
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
inline unsigned int IntersectionAux<dimA,dimB>::component_idx() const
{   return component_idx_; }

template<unsigned int dimA, unsigned int dimB>
inline unsigned int IntersectionAux<dimA,dimB>::is_pathologic() const
{   return pathologic_; }


} // END NAMESPACE
#endif /* INTERSECTIONAUX_H_ */
