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

#include "fem/global_enrichment_func.hh"

#include <armadillo>

template<unsigned int dim, typename S> class SingularityCRTP;

/**
 * CRTP - Curiously recurring template pattern.
 * - avoids virtual member calls for @p val, @p grad, @p div
 */
template<unsigned int dim, template<unsigned int> class S, unsigned int spacedim>
class SingularityCRTP< dim, S<spacedim> > : public GlobalEnrichmentFunc<dim,spacedim>
{
public:
    typedef typename Space<spacedim>::Point Point;
    
    SingularityCRTP(){}
    virtual ~SingularityCRTP(){}
    
    /// scalar enrichment function
    virtual double value(const Point &x) const = 0;
    /// gradient of scalar enrichment function
    virtual Point grad(const Point &x) const = 0;
    
    /// vector enrichment function
    virtual Point vector(const Point &x) const = 0;
    /// divergence of vector enrichment function
    virtual double div(const Point &x) const = 0;
};


///@brief Auxilliary class with geometric data and operations for singularity.
class CircleEllipseProjection
{
public:
    typedef Space<3>::Point Point;
    CircleEllipseProjection(const Point& center, const double& radius,
                            const Point& direction_vector,
                            const Point& normal_vector);
    
    Point center() const
    { return center_;}
    
    double radius() const
    { return radius_;}
    
    Point direction_vector() const
    { return direction_vector_;}
    
    Point normal_vector() const
    { return normal_vector_;}
    
    double circle_area() const
    { return M_PI*radius_*radius_;}
    
    double ellipse_area() const
    { return std::abs(M_PI*a_*b_);}
    
    Point ellipse_a() const
    { return a_*ea_;}
    
    Point ellipse_b() const
    { return b_*eb_;}
    
    /// Projects points from circle plane to ellipse plane.
    void project_to_ellipse_plane(std::vector<Point>& points) const;
    /// Projects points from ellipse plane tot circle plane.
    void project_to_circle_plane(std::vector<Point>& points) const;
    
    /// Returns true if the point @p p lies in the circle (suppose @p p is lying in the circle plane).
    bool point_in_circle(const Point& p) const;
    /// Returns true if the point @p p lies in the ellipse (suppose @p p is lying in the ellipse plane).
    bool point_in_ellipse(const Point& p) const;
    
protected:
    
    Point center_,          ///< Center of the circle.
        direction_vector_,  ///< Direction vector of the singularity (1d element).
        normal_vector_;     ///< Normal vector of the plain (2d element).
        
    double radius_;         ///< Radius of the circle.
    double a_,b_;           ///< Axes sizes of the ellipse.
    Point ea_,eb_;          ///< Unit vectors defining axes of the ellipse.
    
    /// Auxilliary precomputed values for projection; angle between direction and normal vector.
    double cos_a, sin_a, tan_a;
};


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
template<unsigned int spacedim>
class Singularity0D : public SingularityCRTP<2,Singularity0D<spacedim>>
{
public:
    typedef typename Space<spacedim>::Point Point;
    
    Singularity0D(const Point& center, double radius);
    Singularity0D(const Point& center, double radius, const Point& direction_vector, const Point& normal_vector);
    
    //(CRTP functions)
    double value(const Point &x) const override;
    Point grad(const Point &x) const override;
    Point vector(const Point &x) const override;
    double div(const Point &x) const override;

    Point center() const;
    double radius() const;
       
    const CircleEllipseProjection & geometry() const;
    
    /// Computes circumference along edge.
    double circumference() const;
    
    /// Computes @p n points equally distributed along the eqge in 2D ambient space.
    void evaluate_q_points(unsigned int count);
    
    /// Quadrature points on the edge
    const std::vector<Point>& q_points() const;
    
private:
    
    /// center of the singularity
    Point center_;
    
    /// radius of the singularity
    double radius_;
    
    /// points placed around the circle by @p evaluate_q_points function
    std::vector<Point> q_points_;
    
    /// Geometry fucnctionality - plane projections.
    CircleEllipseProjection geom_;
};


template<>
inline Singularity0D<2>::Singularity0D(const Point& center, double radius)
:   center_(center), radius_(radius),
    geom_({center[0],center[1],0},radius, {0,0,1}, {0,0,1})
{
}

template<>
inline Singularity0D<3>::Singularity0D(const Point& center, double radius,
                                       const Point& direction_vector, const Point& normal_vector)
:   center_(center), radius_(radius),
    geom_(center,radius, direction_vector, normal_vector)
{}

template<unsigned int spacedim>
inline typename Singularity0D<spacedim>::Point Singularity0D<spacedim>::center() const
{   return center_;}

template<unsigned int spacedim>
inline double Singularity0D<spacedim>::radius() const
{   return radius_;}


template<unsigned int spacedim>
inline const CircleEllipseProjection & Singularity0D<spacedim>::geometry() const
{   return geom_;}

template<unsigned int spacedim>
inline const std::vector<typename Singularity0D<spacedim>::Point>& Singularity0D<spacedim>::q_points() const
{   return q_points_;}


template<unsigned int spacedim>
double Singularity0D<spacedim>::value(const Point& x) const
{
    double distance = arma::norm(center_-x,2);
    if (distance >= radius_)
        return std::log(distance);
  
    return std::log(radius_);
}

template<unsigned int spacedim>
typename Singularity0D<spacedim>::Point Singularity0D<spacedim>::grad(const Point &x) const
{
    Point grad; //initialize all entries with zero
    grad.zeros();
    
    double distance = arma::norm(center_-x,2);
    if (distance >= radius_)
    {   
        distance = distance * distance;
        grad = (x - center_) / distance;
    }
    return grad;  //returns zero if  (distance <= radius)
}

template<unsigned int spacedim>
typename Singularity0D<spacedim>::Point Singularity0D<spacedim>::vector(const Point &x) const
{
    return this->grad(x);
}

template<unsigned int spacedim>
double Singularity0D<spacedim>::div(const Point& x) const
{
    return 0;
}

template<unsigned int spacedim>
inline double Singularity0D<spacedim>::circumference() const
{
    return 2 * M_PI * radius_;
}

template<>
inline void Singularity0D<2>::evaluate_q_points(unsigned int count)
{   
    //q_points_.clear();
    q_points_.resize(count);
    
    double phi = 2*M_PI / count;
    double offset = 1e-10;
    
    for(unsigned int i=0; i < count; i++)
        q_points_[i] = {center_[0] + radius_*std::cos(i*phi+offset),
                        center_[1] + radius_*std::sin(i*phi+offset)};
}

template<>
inline void Singularity0D<3>::evaluate_q_points(unsigned int count)
{
    //q_points_.clear();
    q_points_.resize(count);  

    double phi = 2*M_PI / count;
    double offset = 1e-10;
    
    for(unsigned int i=0; i < count; i++)
        q_points_[i] = center_ + std::cos(i*phi+offset)*geom_.ellipse_b() + std::sin(i*phi+offset)*geom_.ellipse_a();
}


#endif // SINGULARITY_HH_