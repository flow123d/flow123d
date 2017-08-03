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

// template<unsigned int dim, class S> class SingularityCRTP;

/**
 * CRTP - Curiously recurring template pattern.
 * - avoids virtual member calls for @p val, @p grad, @p div
 */
template<unsigned int dim, class S>
class SingularityCRTP : public GlobalEnrichmentFunc<dim,3>
{
public:
    typedef Space<3>::Point Point;
    
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
    
    double distance(const Point &p) const
    { return arma::norm(center_-p,2);}
    
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
class Singularity0D : public SingularityCRTP<2,Singularity0D>
{
public:
    typedef Space<3>::Point Point;
    
//     Singularity0D(const Point& center, double radius);
    Singularity0D(const Point& center, double radius, const Point& direction_vector, const Point& normal_vector);
    
    //(CRTP functions)
    double value(const Point &x) const override;
    Point grad(const Point &x) const override;
    Point vector(const Point &x) const override;
    double div(const Point &x) const override;

    Point center() const
    { return geom_.center();}
    double radius() const
    { return geom_.radius();}
       
    const CircleEllipseProjection & geometry() const
    { return geom_;}
    
    /// Computes circumference along edge.
    double circumference() const;
    
    /// Computes @p n points equally distributed along the eqge in 2D ambient space.
    void evaluate_q_points(unsigned int count);
    
    /// Quadrature points on the edge
    const std::vector<Point>& q_points() const
    { return q_points_;}
    
    /// Value of the permeability coefficient between dimensions.
    double sigma() const
    { return sigma_;}
    void set_sigma(double s)
    { sigma_ = s;}
    
    double pressure() const
    { return pressure_;}
    void set_pressure(double p)
    { pressure_ = p;}
    
private:
    
    /// points placed around the circle by @p evaluate_q_points function
    std::vector<Point> q_points_;
    
    /// Geometry fucnctionality - plane projections.
    CircleEllipseProjection geom_;
    
    /// Value of the permeability coefficient between dimensions.
    double sigma_;
    double pressure_;
    
    double radius_rounding_low_bound_;
};


// inline Singularity0D::Singularity0D(const Point& center, double radius)
// :   center_(center), radius_(radius),
//     geom_({center[0],center[1],0},radius, {0,0,1}, {0,0,1})
// {
//     radius_rounding_low_bound_ = radius - 1e-8*radius;
// }

inline double Singularity0D::value(const Point& x) const
{
    double distance = geom_.distance(x);
    if (distance >= radius_rounding_low_bound_)
        return std::log(distance);
  
    return std::log(geom_.radius());
}

inline Singularity0D::Point Singularity0D::grad(const Point &x) const
{
    Point grad; //initialize all entries with zero
    grad.zeros();
    
    double distance = geom_.distance(x);
    if (distance >= radius_rounding_low_bound_)
    {   
        distance = distance * distance;
        grad = (x - geom_.center()) / distance;
    }
    return grad;  //returns zero if  (distance <= radius)
}

inline Singularity0D::Point Singularity0D::vector(const Point &x) const
{
    static const double t = -0.5/M_PI;
    return t * this->grad(x);
}

inline double Singularity0D::div(const Point& x) const
{
    return 0;
}

inline double Singularity0D::circumference() const
{
    //TODO: https://en.wikipedia.org/wiki/Ellipse#Circumference
    return 2 * M_PI * geom_.radius();
}






























///@brief Auxilliary class with geometric data and operations for singularity.
class CylinderGeometry
{
public:
    typedef Space<3>::Point Point;
    CylinderGeometry(const Point& a, const Point& b,
                     const double& radius);
    
    Point a() const
    { return a_;}
    Point b() const
    { return b_;}
    
    double radius() const
    { return radius_;}
    
    Point direction_vector() const
    { return direction_vector_;}
     
    /// Distance vector from point p to the axis.
    Point dist_vector(const Point &p) const;
    
    double distance(const Point &p) const
    { return arma::norm(dist_vector(p),2);}
    
    /// Returns true if the point @p p lies in the circle (suppose @p p is lying in the circle plane).
    bool point_in_cylinder(const Point& p) const
    { return distance(p) <= radius_;}
    
    /// Computes lateral surface of the cylinder.
    double lateral_surface() const
    { return 2*M_PI*radius_*arma::norm(direction_vector_,2);}
    
    /// Computes volume of the cylinder.
    double volume() const
    { return M_PI * radius_ * radius_ * arma::norm(direction_vector_,2);}
    
protected:
    
    Point a_,b_,            ///< End points of the axis of the cylinder.
        direction_vector_;  ///< Direction vector of the axis.
        
    double radius_;         ///< Radius of the cylinder base.
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
class Singularity1D : public SingularityCRTP<3,Singularity1D>
{
public:
    typedef Space<3>::Point Point;
    
    Singularity1D(const Point& a, const Point& b, double radius);
    
    //(CRTP functions)
    double value(const Point &x) const override;
    Point grad(const Point &x) const override;
    Point vector(const Point &x) const override;
    double div(const Point &x) const override;
    
    /// Returns geometry object.
    const CylinderGeometry & geometry() const
    { return geom_;}
    
    /** Computes @p n*m points equally distributed on the lateral surface of cylinder.
     * @param n is number of points around the cylinder in one layer
     * @param m is number of layers of points in z coordinate of the cylinder
     */
    void evaluate_q_points(unsigned int n, unsigned int m);
    
    /// Quadrature points on the edge
    const std::vector<Point>& q_points() const
    { return q_points_;}
    
    /// Value of the permeability coefficient between dimensions.
    double sigma() const
    { return sigma_;}
    void set_sigma(double s)
    { sigma_ = s;}
    
    double pressure() const
    { return pressure_;}
    void set_pressure(double p)
    { pressure_ = p;}
    
private:
    /// points placed around the circle by @p evaluate_q_points function
    std::vector<Point> q_points_;
    
    /// Geometry fucnctionality - plane projections.
    CylinderGeometry geom_;
    
    /// Value of the permeability coefficient between dimensions.
    double sigma_;
    double pressure_;
    
    double radius_rounding_low_bound_;
};


inline double Singularity1D::value(const Point& x) const
{
    double distance = geom_.distance(x);
    if (distance >= radius_rounding_low_bound_)
        return std::log(distance);
  
    return std::log(geom_.radius());
}

inline Space<3>::Point Singularity1D::grad(const Point &x) const
{
    Point grad = geom_.dist_vector(x);    
    double distance = arma::norm(grad,2);
    
    if (distance >= radius_rounding_low_bound_)
    {   
        distance = distance * distance;
        grad = (grad) / distance;
    }
    else{
        grad.zeros();
    }
    return grad;  //returns zero if  (distance <= radius)
}

inline Space<3>::Point Singularity1D::vector(const Point &x) const
{
    static const double t = -0.5/M_PI;
    return t * this->grad(x);
}

inline double Singularity1D::div(const Point& x) const
{
    return 0;
}

#endif // SINGULARITY_HH_