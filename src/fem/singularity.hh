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
 * @file    singularity.hh
 * @brief   Definition of singularities for enrichment in XFEM.
 * @author  Pavel Exner
 */

#ifndef SINGULARITY_HH_
#define SINGULARITY_HH_

#include "system/global_defs.h"

#include "mesh/point.hh"
#include "mesh/elements.h"

#include <armadillo>

/** @brief Point singularity in 2D plane in 2D or 3D ambient space.
 * 
 * Singularity is defined by its center point and radius.
 * Radius is the actual dimension of the singularity.
 * 
 * It provides global enrichment function @p val and its gradient @p grad for XFEM.
 * 
 * It can provide points on the circumference of the singularity edge.
 * These must be explicitely computed by @p evaluate_q_points and later get by @p q_points.
 * 
 * In 3D ambient space: 
 * - the plane in which the singularity lies must be set for evaluation of the @p q_points
 * - the plane is defined by the element in which the singularity lies.
 * - for simplicity, we do not map circle to ellipse in the plane
 *   (we suppose always circular cross section of 1D-2D)
 * 
 */
template<int spacedim>
class Singularity0D
{
public:
    typedef typename Space<spacedim>::Point Point;
    
    Singularity0D(const Point& center, double radius);
    
    Point center() const;
    double radius() const;
    
    double val(const Point &x) const;
    Point grad(const Point &x) const;
    
    /// Computes circumference along edge.
    double circumference() const;
    
    /// Computes @p n points equally distributed along the eqge in 2D ambient space.
    void evaluate_q_points(unsigned int count);
    
    /** Computes @p n points equally distributed along the eqge in 3D ambient space.
     * 
     * @param ele - element in which the singularity lies; defines the plane in 3D ambient space.
     */
    void evaluate_q_points(unsigned int count, ElementFullIter ele);
    
    /// Quadrature points on the edge
    const std::vector<Point>& q_points() const;
    
private:
    /// center of the singularity
    Point center_;
    
    /// radius of singularity
    double radius_;
    
    /// points placed around the circle by @p evaluate_q_points function
    std::vector<Point> q_points_;
};

template<int spacedim>
Singularity0D<spacedim>::Singularity0D(const Point& center, double radius)
:   center_(center), radius_(radius)
{}

template<int spacedim>
inline typename Singularity0D<spacedim>::Point Singularity0D<spacedim>::center() const
{   return center_;}

template<int spacedim>
inline double Singularity0D<spacedim>::radius() const
{   return radius_;}

template<int spacedim>
inline const std::vector<typename Singularity0D<spacedim>::Point>& Singularity0D<spacedim>::q_points() const
{   return q_points_;}


template<int spacedim>
double Singularity0D<spacedim>::val(const Point& x) const
{
    double distance = arma::norm(center_-x,2);
    if (distance >= radius_)
        return std::log(distance);
  
    return std::log(radius_);
}


template<int spacedim>
typename Singularity0D<spacedim>::Point Singularity0D<spacedim>::grad(const Point &x) const
{ 
    Point grad; //initialize all entries with zero
    grad.zeros();
    
    double distance = arma::norm(center_-x,2);
    if (distance >= radius_)
    {   
        distance = distance * distance;
        grad = (x-center_) / distance;
    }
    return grad;  //returns zero if  (distance <= radius)
}

template<int spacedim>
inline double Singularity0D<spacedim>::circumference() const
{
    return 2*M_PI*radius_;
}

template<>
inline void Singularity0D<2>::evaluate_q_points(unsigned int count)
{   
    //q_points_.clear();
    q_points_.resize(count);
    
    double phi = 2*M_PI / count;
    double offset = 1e-10;
    
    for(unsigned int i=0; i < count; i++)
        q_points_[i] = {center_[0]+radius_*std::cos(i*phi+offset),
                        center_[1]+radius_*std::sin(i*phi+offset)};
}

template<>
inline void Singularity0D<3>::evaluate_q_points(unsigned int count, ElementFullIter ele)
{
    ASSERT_DBG(ele->dim() == 2);
    
    //q_points_.clear();
    q_points_.resize(count);
 
    Point u;   //first unit directional vector in the plane of singularity
    Point v;   //second unit directional vector in the plane of singularity
    Point n;   //normal to element
    //directional vectors of triangle sides
    Point e1 = ele->node[1]->point() - ele->node[0]->point();
    Point e2 = ele->node[2]->point() - ele->node[0]->point();
    n = arma::cross(e1,e2);
    v = arma::cross(n,e1);
    u = e1/arma::norm(e1,2);
    v = v/arma::norm(v,2);
    
    double phi = 2*M_PI / count;
    double offset = 1e-10;
    
    for(unsigned int i=0; i < count; i++)
        q_points_[i] = center_ + radius_*(std::cos(i*phi+offset)*u + std::sin(i*phi+offset)*v);
}


#endif // SINGULARITY_HH_