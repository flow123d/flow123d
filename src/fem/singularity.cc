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
 * @file    singularity.cc
 * @brief   Definition of geometry functionalities of singularities.
 * @author  Pavel Exner
 */

#include "singularity.hh"
#include "system/global_defs.h"

#include <armadillo>


CircleEllipseProjection::CircleEllipseProjection(const CircleEllipseProjection::Point& center,
                                                 const double& radius,
                                                 const CircleEllipseProjection::Point& direction_vector,
                                                 const CircleEllipseProjection::Point& normal_vector)
: center_(center), radius_(radius)
{
    //normalize vectors
    direction_vector_ = direction_vector / arma::norm(direction_vector,2);
    normal_vector_ = normal_vector / arma::norm(normal_vector,2);
    Point& u = direction_vector_;
    Point& n = normal_vector_;
    
    //cosine of the angle between n and u
    cos_a = arma::dot(n,u);
    ASSERT(abs(cos_a) > std::numeric_limits<double>::epsilon());
    
    sin_a = std::sqrt(1-cos_a*cos_a);
    tan_a = sin_a / cos_a;
    
    // ellipse axes sizes
    a_ = radius;
    b_ = radius / cos_a;
    
    // ellipse axes unit vectors
    if(sin_a > std::numeric_limits<double>::epsilon()){
        ea_ = arma::cross(n,u);
        eb_ = arma::cross(n,ea_);
    }
    else { 
        //direction and normal vectors are parallel
        //we make the circle axes using the z axis
        Point z = {0,0,1};
        double cos_z = arma::dot(z,u);
        if( abs(1-cos_z) > std::numeric_limits<double>::epsilon()){
            ea_ = arma::cross(u,z);
            eb_ = arma::cross(u,ea_);
        }
        else {
            //direction vector is parallel to z axis => use x and y axes
            ea_ = {1,0,0};
            eb_ = {0,1,0};
        }
    }
    
    //normalize the vectors
    ea_ = ea_ / arma::norm(ea_,2);
    eb_ = eb_ / arma::norm(eb_,2);
    
//     (center + u).print(DebugOut(),"direction");
//     (center + n/arma::norm(n,2)).print(DebugOut(),"normal");
//     center.print(DebugOut(),"center");
//     (center+a_*ea_).print(DebugOut(),"a");
//     (center+b_*eb_).print(DebugOut(),"b");
}

bool CircleEllipseProjection::point_in_circle(const CircleEllipseProjection::Point& p) const
{
    // distance of point is less than radius
    return arma::norm(center_ - p,2) <= radius_ + 8*std::numeric_limits<double>::epsilon();
}

bool CircleEllipseProjection::point_in_ellipse(const CircleEllipseProjection::Point& p) const
{
    // Point of ellipse is defined with angle parameter t: 
    // p = c + a * cost * ea + b * sint * eb;
    // where ea and eb are unit vectors defining axes of ellipse
    // denoting: 
    Point rhs = p - center_;
    Point u = a_ * ea_,
          v = b_ * eb_;
    // we make 2 equations for c (cost) and s (sint):
    // rhs.u = c*u.u + s*v.u;
    // rhs.v = c*u.v + s*v.v;
    
    // c and s are coordinates in the u,v coordinate system of the ellipse:
    // (c*a)^2 / a^2 + (s*b)^2 / b^2 = 1
    // and we want the point inside, therefore:
    // c^2 + s^2 <= 1
    
    // dot products:
    double
        d00 = arma::dot(u,u),
        d01 = arma::dot(u,v),
        d02 = arma::dot(u,rhs),
        d11 = arma::dot(v,v),
        d12 = arma::dot(v,rhs);
    double invdenom = 1.0/ (d00*d11 - d01*d01);
    double c = (d02*d11 - d12*d01) * invdenom,
           s = (d12*d00 - d02*d01) * invdenom;
    
    return c*c + s*s <= 1 + 8*std::numeric_limits<double>::epsilon();
}


void CircleEllipseProjection::project_to_ellipse_plane(vector< CircleEllipseProjection::Point >& points) const
{
    //direction and normal vectors are parallel => perfect circle
    if(sin_a < std::numeric_limits<double>::epsilon()) return;
    
    // ww is a vector in the plane of circle - we later compute angle in the circle plane from ww
    Point ww = arma::cross(direction_vector_,ea_);
    ww = ww / arma::norm(ww,2);
    
    for(unsigned int j=0; j < points.size(); j++){
        Point r = points[j] - center_;
        double r_norm = arma::norm(r,2);
        double cos_r = arma::dot(r,ww) / r_norm;
    
        //maximal translation is in ww direction
        //the factor of translation is given by cos of angle between r and ww
        double x = tan_a * r_norm;
        points[j] += - cos_r * x * direction_vector_;
    }
}


void CircleEllipseProjection::project_to_circle_plane(vector< CircleEllipseProjection::Point >& points) const
{
    //direction and normal vectors are parallel => perfect circle
    if(sin_a < std::numeric_limits<double>::epsilon()) return;
    
    for(unsigned int j=0; j < points.size(); j++){
        Point r = points[j] - center_;
        double r_norm = arma::norm(r,2);
        double cos_r = arma::dot(r,eb_) / r_norm;
       
        //maximal translation is in eb direction
        //the factor of translation is given by cos of angle between r and eb
        double x = sin_a * r_norm;
        points[j] += cos_r * x * direction_vector_;
    }
}