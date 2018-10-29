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
 * @file    global_enrichment_function.hh
 * @brief   Definition of the interface for global enrichment function in XFEM.
 * @author  Pavel Exner
 */

#ifndef GLOBAL_ENRICHMENT_FUNCTION_HH_
#define GLOBAL_ENRICHMENT_FUNCTION_HH_

#include "mesh/point.hh"

class Geometry;

/** @brief Abstract class defining the interface of global enrichment function.
 * 
 * Enrichment done by local solution must provide @p val (and @p grad if needed) functions for scalar enrichment
 * and @p div_grad (and @p grad_grad if needed) for vectorial enrichment.
 */
template<int dim, int spacedim>
class GlobalEnrichmentFunc
{
public:
    typedef typename Space<spacedim>::Point Point;
    
    GlobalEnrichmentFunc(){}
    virtual ~GlobalEnrichmentFunc(){}
    
    /// scalar enrichment function
    virtual double value(const Point &x) const = 0;
    /// gradient of scalar enrichment function
    virtual Point grad(const Point &x) const = 0;
    
    /// vector enrichment function
    virtual Point vector(const Point &x) const = 0;
//     virtual arma::mat::fixed<spacedim,spacedim> grad_grad(const Point &x) const = 0;
    /// divergence of vector enrichment function
    virtual double div(const Point &x) const = 0;
    
    virtual const Geometry& geometry() const = 0;
};

class Geometry
{
    typedef typename Space<3>::Point Point;
public:
    virtual Point dist_vector(const Point &p) const = 0;
    virtual double distance(const Point &p) const = 0;
    virtual double effective_surface() const = 0;
    virtual double volume() const = 0;
    virtual bool point_inside(const Point &p) const = 0;
    const std::vector<Point>& q_points() const
    { return q_points_;};
    
protected:
    /// quadrature points on effective surface computed by @p evaluate_q_points function
    std::vector<Point> q_points_;
};

#endif // GLOBAL_ENRICHMENT_FUNCTION_HH_