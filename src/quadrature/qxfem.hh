/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    qxfem.hh
 * @brief   Adaptive quadrature rules for usage in XFEM.
 * @author  Pavel Exner
 */

#ifndef QXFEM_HH_
#define QXFEM_HH_

#include "system/global_defs.h"

#include "quadrature/quadrature.hh"

#include "mesh/point.hh"

template<int dim, int spacedim> class QXFEMFactory;

/** @brief Class representing quadrature for XFEM, with adaptively distributed points.
 * Includes quadrature points both in real coordinates and on the reference element.
 */
template<int dim, int spacedim>
class QXFEM : public Quadrature<dim> {
public:
    typedef typename Space<spacedim>::Point Point;
    /// Empty constructor
    QXFEM(){}
    
    /// @name Getters.
    //@{
    
    /// Getter for vector of quadrature points in real coordinates.
    const std::vector<Point> & get_real_points() const;
    
    /// Getter for i-th quadrature point in real coordinates.
    const Point & real_point(unsigned int i) const;
    //@}
    
private:
    /// Vector of quadrature points in real coordinates.
    std::vector<Point> real_points_;
    
    friend class QXFEMFactory<dim,spacedim>;
};


template<int dim, int spacedim>
inline const typename Space<spacedim>::Point& QXFEM<dim, spacedim>::real_point(unsigned int i) const
{   ASSERT_DBG(i < real_points_.size());
    return real_points_[i];
}

template<int dim, int spacedim>
inline const std::vector< typename Space<spacedim>::Point >& QXFEM<dim, spacedim>::get_real_points() const
{   return real_points_;}



#endif // QXFEM_HH_