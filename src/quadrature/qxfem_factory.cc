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
 * @file    qxfem_factory.cc
 * @brief   XFEM adaptive quadrature factory.
 * @author  Pavel Exner
 */

#include "system/global_defs.h"

#include "qxfem_factory.hh"

#include "quadrature/quadrature_lib.hh"

#include "mesh/ref_element.hh"
#include "mesh/elements.h"

#include "fem/mapping_p1.hh"

const double QXFEMFactory::distance_criteria_factor_2d_ = 0.25;
const double QXFEMFactory::distance_criteria_factor_3d_ = 1.0;

void QXFEMFactory::clear()
{
    simplices_.clear();
    level_ = 0;
    level_offset_ = 0;
}

template<int dim>
double QXFEMFactory::AuxSimplex::compute_max_h() const
{
    ASSERT(nodes.size() == dim+1);
    double max_h = 0;
    for(unsigned int li=0; li<RefElement<dim>::n_lines; li++){
        Point v = nodes[RefElement<dim>::interact(Interaction<0,1>(li))[0]]
                 -nodes[RefElement<dim>::interact(Interaction<0,1>(li))[1]];
        double t = arma::norm(v,2);
        if(t>max_h) max_h = t;
    }
    return max_h;
}

std::shared_ptr< QXFEM< 2,3 > > QXFEMFactory::create_singular(
                                                const std::vector<Singularity0DPtr> & sing,
                                                ElementFullIter ele)
{
    clear();
    std::shared_ptr<QXFEM<2,3>> qxfem = make_shared<QXFEM<2,3>>();
    
    //create first auxsimplex:
    AuxSimplex s;
    for(unsigned int i=0; i < ele->n_nodes(); i++)
        s.nodes.push_back(ele->node[i]->point());
    
    //we suppose counterclockwise node order
    Point v0 = s.nodes[1] - s.nodes[0],   // 0. edge of triangle
          v1 = s.nodes[2] - s.nodes[1];   // 1. edge of triangle
    if(v0[1]*v1[2] - v0[2]*v1[1] + v0[2]*v1[0] - v0[0]*v1[2] + v0[0]*v1[1] - v0[1]*v1[0] < 0)
    {
        DebugOut() << "change node order\n";
        std::swap(s.nodes[0],s.nodes[1]);
    }
    
    //TODO
    // WORKS only for parallel wells
    // ideas:
    // 1) project to circle plane -> refine -> project to another circle plane -> refine -> ...
    // 2) work with ellipses
    
    // project the simplex into the plane of the circle
    sing[0]->geometry().project_to_circle_plane(s.nodes);
    
    max_h_level_0_ = s.compute_max_h<2>();
    s.active = true;
    simplices_.push_back(s);
    
    for(unsigned int i=0; i < max_level_; i++)
    {
        level_ = i;
        DebugOut() << "QXFEM level: " << level_ << "\n";
        unsigned int n_to_refine = mark_refinement_2D(sing);
        if(n_to_refine == 0) break; // end refinement
        
        refine_level<2>(n_to_refine);
    }
    mark_refinement_2D(sing);
    
    // project the simplices back to the plane of the ellipse
    for(AuxSimplex& simp : simplices_){
        sing[0]->geometry().project_to_ellipse_plane(simp.nodes);
    }
    
    distribute_qpoints<2>(qxfem->real_points_, qxfem->weights, sing, ele->measure());
    map_real_to_unit_points<2>(qxfem->real_points_, qxfem->quadrature_points, ele);
    DebugOut().fmt("n_real_qpoints {} {}\n",qxfem->real_points_.size(), qxfem->quadrature_points.size());
    
    return qxfem;
}


std::shared_ptr< QXFEM< 3,3 > > QXFEMFactory::create_singular(
                                                const std::vector<Singularity1DPtr> & sing,
                                                ElementFullIter ele)
{
    clear();
    std::shared_ptr<QXFEM<3,3>> qxfem = make_shared<QXFEM<3,3>>();
    
    //create first auxsimplex:
    AuxSimplex s;
    for(unsigned int i=0; i < ele->n_nodes(); i++)
        s.nodes.push_back(ele->node[i]->point());
    
//     //TODO: revise for tetrahedron
//     //we suppose counterclockwise node order
//     Point v0 = s.nodes[1] - s.nodes[0],   // 0. edge of triangle
//           v1 = s.nodes[2] - s.nodes[1];   // 1. edge of triangle
//     if(v0[1]*v1[2] - v0[2]*v1[1] + v0[2]*v1[0] - v0[0]*v1[2] + v0[0]*v1[1] - v0[1]*v1[0] < 0)
//     {
//         DebugOut() << "change node order\n";
//         std::swap(s.nodes[0],s.nodes[1]);
//     }
    
    max_h_level_0_ = s.compute_max_h<3>();
    s.active = true;
    simplices_.push_back(s);
    
    for(unsigned int i=0; i < max_level_; i++)
    {
        level_ = i;
        DebugOut() << "QXFEM level: " << level_ << "\n";
        unsigned int n_to_refine = mark_refinement_3D(sing);
        if(n_to_refine == 0) break; // end refinement
        
        refine_level<3>(n_to_refine);
    }
    mark_refinement_3D(sing);
    
    distribute_qpoints<3>(qxfem->real_points_, qxfem->weights, sing, ele->measure());
    map_real_to_unit_points<3>(qxfem->real_points_, qxfem->quadrature_points, ele);
    DebugOut().fmt("n_real_qpoints {} {}\n",qxfem->real_points_.size(), qxfem->quadrature_points.size());
    
    return qxfem;
}


unsigned int QXFEMFactory::mark_refinement_2D(const std::vector<Singularity0DPtr> & sing)
{    
    unsigned int n_simplices_to_refine = 0; 
  /* there are several cases that can happen:
   * 1] a node of a square can be inside a well       -> refine
   * 2] if all nodes of a square are inside a well 
   *    the whole square is inside                    -> no refine
   * 3] if no node of a square is inside a well 
   *    and the center of a well is inside the square 
   *    then the whole well is inside the square      -> refine
   * 4] if no nodes are inside are in a well
   *    and the well edge crosses the square line     -> refine
   * 5] if the minimum distance of square vertex 
   *    to well edge is smaller than square diameter  -> refine
   */
  
    // go through simplices on the current level
    for(unsigned int i = level_offset_; i < simplices_.size(); i++)
    {
        AuxSimplex& s = simplices_[i];
        for(unsigned int j = 0; j < sing.size(); j++)
        {
//             DebugOut().fmt("QXFEM test simplex {}, singularity {}\n",i,j);
            double distance_sqr = -1;
            double max_h = max_h_level_0_/(1<<level_);
            int res = sigularity0D_distance(*sing[j],s, distance_sqr);
            if(res > 0) {
                s.refine = true;
                s.sing_id = j;
                n_simplices_to_refine++;
                break;
            }
            // distance criterion
            if(distance_sqr > 0)
            {
                double rmin = std::sqrt(distance_sqr)-sing[j]->radius();
                rmin = rmin*rmin;
                if( max_h > distance_criteria_factor_2d_ * rmin)
                {
                    s.refine = true;
                    n_simplices_to_refine++;
                    break;
                }
            }
            if(res == 0) {
                s.active = false;
                break;
            }
        }
    }
    
    return n_simplices_to_refine;
}


unsigned int QXFEMFactory::mark_refinement_3D(const std::vector<Singularity1DPtr> & sing)
{    
    unsigned int n_simplices_to_refine = 0; 
  
    // go through simplices on the current level
    for(unsigned int i = level_offset_; i < simplices_.size(); i++)
    {
        AuxSimplex& s = simplices_[i];
        for(unsigned int j = 0; j < sing.size(); j++)
        {
//             DebugOut().fmt("QXFEM test simplex {}, singularity {}\n",i,j);
            double distance = -1;
            double max_h = max_h_level_0_/(1<<level_);
            int res = singularity1D_distance(*sing[j],s, distance);
//             DBGVAR(res);
            if(res > 0) {
                s.refine = true;
                s.sing_id = j;
                n_simplices_to_refine++;
                break;
            }
            // distance criterion
            if(distance > 0)
            {
                double rmin = distance-sing[j]->geometry().radius();
                rmin = rmin*rmin;
                rmin = rmin;
                if( max_h > distance_criteria_factor_3d_ * rmin)
                {
                    s.refine = true;
                    n_simplices_to_refine++;
                    break;
                }
            }
            if(res == 0) {
                s.active = false;
                break;
            }
        }
    }
    
    return n_simplices_to_refine;
}

int QXFEMFactory::sigularity0D_distance(const Singularity0D& w,
                                        const AuxSimplex& s,
                                        double& distance_sqr)
{
    ASSERT_DBG(s.nodes.size() == 3);
    
    //we suppose counterclockwise node order
    
//             DebugOut() << "QXFEM refine test1\n";
            // TEST 1: Vertex within circle
            // http://www.phatcode.net/articles.php?id=459
            Point c0 = w.center() - s.nodes[0],
                  c1 = w.center() - s.nodes[1],
                  c2 = w.center() - s.nodes[2];
            
            
            double radius_sqr = w.radius() * w.radius();
            double c0sqr = arma::dot(c0,c0) - radius_sqr;
            double c1sqr = arma::dot(c1,c1) - radius_sqr;
            double c2sqr = arma::dot(c2,c2) - radius_sqr;
            
            int count_nodes = 0;
            if (c0sqr <= 0) count_nodes++;
            if (c1sqr <= 0) count_nodes++;
            if (c2sqr <= 0) count_nodes++;
            
//             DebugOut() << "count_nodes " << count_nodes << "\n";
            if(count_nodes > 0 ){
                if(count_nodes == 3) { //triangle inside circle
                    return 0;
                }
                else {
                    return 1;
                }
            }
            
//             DebugOut() << "QXFEM refine test2\n";
            // TEST 2: Circle centre within triangle
            // http://www.blackpawn.com/texts/pointinpoly/
            // using barycentric coordinates test
            Point v0 = s.nodes[1] - s.nodes[0],   // 0. edge of triangle
                  v1 = s.nodes[2] - s.nodes[0];   // 1. edge of triangle

            // v2 = c0
            // Point in triangle is defined: 
            // v2 = s*v0 + t*v1;
            // make 2 equations for s and t:
            // v2.v0 = s*v0.v0 + t*v1.v0;
            // v2.v1 = + s*v0.v1 + t*v1.v1;
            
            // dot products:
            double
                d00 = arma::dot(v0,v0),
                d01 = arma::dot(v0,v1),
                d02 = arma::dot(v0,c0),
                d11 = arma::dot(v1,v1),
                d12 = arma::dot(v1,c0);
            double invdenom = 1.0/ (d00*d11 - d01*d01);
            double t0 = (d02*d11 - d12*d01) * invdenom,
                   t1 = (d12*d00 - d02*d01) * invdenom;
            
//             DebugOut().fmt("t0 = {} t1 = {}\n", t0, t1);
            
            if ( (t0 >= 0) && (t1 >= 0) && (t0+t1 <= 1)){
                return 2;
            }
            
//             DebugOut() << "QXFEM refine test3\n";
            // ; TEST 3: Circle intersects edge
            
            double k0 = arma::dot(c0,v0);

            if( (k0 >= 0) && (k0 <= d00) && (c0sqr * d00 <= k0*k0)){
                return 3;
            }

            double k1 = arma::dot(c0,v1);
            if( (k1 >= 0) && (k1 <= d11) && (c0sqr * d11 <= k1*k1)){
                return 4;
            }
            
            Point v2 = s.nodes[2] - s.nodes[1];   // 2. edge of triangle
            double d22 = arma::dot(v2,v2);
            double k2 = arma::dot(c1,v2);
            if( (k2 >= 0) && (k2 <= d22) && (c1sqr * d22 <= k2*k2)){
                return 5;
            }
//             DebugOut() << "QXFEM simplex not refined\n";
            
//             DebugOut().fmt("d00 = {}\n", d00);
//             DebugOut().fmt("d11 = {}\n", d11);
//             DebugOut().fmt("d22 = {}\n", d22);
//             DebugOut().fmt("k0 = {}  {}\n",k0, k0*k0/d00);
//             DebugOut().fmt("k1 = {}  {}\n",k1, k1*k1/d11);
//             DebugOut().fmt("k2 = {}  {}\n",k2, k2*k2/d22);
            
            // possibly compute distance_sqr
            if(distance_sqr < 0) 
            {
                if( (k0 >= 0) && (k0 <= d00) && (t1 < 0)){
//                     DebugOut() << "0 edge\n";
                    distance_sqr = arma::dot(c0,c0) - k0*k0/d00;
                    return -1;
                }

                if( (k1 >= 0) && (k1 <= d11) && (t0 < 0) ){
//                     DebugOut() << "1 edge\n";
                    distance_sqr = arma::dot(c0,c0) - k1*k1/d11;
                    return -1;
                }
                
                if( (k2 >= 0) && (k2 <= d22) && (t0+t1 > 1) ){
//                     DebugOut() << "2 edge\n";
                    distance_sqr = arma::dot(c1,c1) - k2*k2/d22;
                    return -1;
                }
                
                if( (k0 < 0) && (k1 < 0) ){
//                     DebugOut() << "v0\n";
                    distance_sqr = arma::dot(c0,c0);
                    return -1;
                }
                if( (k0 > 0) && (k2 < 0) ){
//                     DebugOut() << "v1\n";
                    distance_sqr = arma::dot(c1,c1);
                    return -1;
                }
                if( (k1 > 0) && (k2 > 0) ){
//                     DebugOut() << "v2\n";
                    distance_sqr = arma::dot(c2,c2);
                    return -1;
                }
            }
                
            return -1;
}


// int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
// {
//   int i, j, c = 0;
//   for (i = 0, j = nvert-1; i < nvert; j = i++) {
//     if ( ((verty[i]>testy) != (verty[j]>testy)) &&
//      (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
//        c = !c;
//   }
//   return c;
// }

int QXFEMFactory::singularity1D_distance(const Singularity1D& w,
                                            const QXFEMFactory::AuxSimplex& s,
                                            double& distance)
{
    const unsigned int n_nodes = s.nodes.size();
    ASSERT_DBG(s.nodes.size() == n_nodes);
    
    //TODO idea: keep results of mutual position for AuxSimplex and on refined levels check only relevant ones..
    const CylinderGeometry& geom = w.geometry();
    
    //TEST 1: nodes of tetrahedron in cylinder

    std::vector<Point> dist_v(n_nodes);
    unsigned int count_nodes = 0;
    double rsqr = geom.radius() * geom.radius();
    for(unsigned int i=0; i<n_nodes; i++){
        dist_v[i] = geom.dist_vector(s.nodes[i]);
        if(arma::dot(dist_v[i],dist_v[i]) < rsqr) count_nodes++;
    }
    
//     DBGVAR(count_nodes);
    if(count_nodes > 0 ){
        distance = 0.0;
        if(count_nodes == n_nodes) return 0; //tetrahedron inside cylinder, do not refine
        else return 1;
    }
    
    //TEST 2: check sing inside tetrahedron
    //TODO: consider line segment, not entire line !!!
    
    Point u = geom.direction_vector();
    Point un = u / arma::norm(u,2);
    
    // projection to plane with un normal vector and A origin
    arma::mat proj;
    proj.eye(3,3);
    proj = proj - un*un.t();
    
    std::vector<Point> v(n_nodes);
    for(unsigned int i=0; i<n_nodes; i++)
        v[i] = proj * (s.nodes[i] - geom.a());
    
    BoundingBox bp (v);
    if(bp.contains_point(Point({0,0,0}))){
        distance = 0.0;
        return 2;
    }
    else{
        Point c = bp.center();
        double max_axis_size = 0;
        for(unsigned int i=0; i<3; i++)
            if(bp.size(i) > max_axis_size) max_axis_size = bp.size(i);
        
        // perpendicular distance from line to BB center
        double d = arma::norm(c,2);
        
        if(d > (geom.radius())) {
            distance = d;
            return -1;
        }
        else{
            distance = 0.0;
            return 3;
        }
    }
    
    
    /* TEST 3: compute distance between S center and line w
     *    if distance >= M*S_max or M*radius, where M is constant (e.g. 3), then use this distance
     *    (or use BoundingBox)
     */
    
//     BoundingBox b(s.nodes);
//     Point c = b.center();
//     
//     double max_axis_size = 0;
//     for(unsigned int i=0; i<3; i++)
//         if(b.size(i) > max_axis_size) max_axis_size = b.size(i);
//     
//     // perpendicular distance from line to BB center
//     double d = geom.distance(c);
//     
//     max_h = max_axis_size;
//     if(d > (geom.radius())) {
//         distance = d;
//         return -1;
//     }
//     else{
//         distance = 0.0;
//         return 4;
//     }
//     if(d > (max_axis_size + geom.radius())) {
//         distance = d;
//         max_h = max_axis_size;
//         return -1;
//     }
//     else {
//         distance = 0.0;
//         return 1;
//     }
    
    
    
    //find minimal axis    
//     BoundingBox b (s.nodes);
//     double min_axis_size = std::numeric_limits<double>::max();
//     unsigned int min_axis=0;
//     for(unsigned int i=0; i<3; i++)
//         if(b.size(i) < min_axis_size){
//             min_axis_size = b.size(i);
//             min_axis = i;
//         }
//     double min_axis_size = 0;
//     unsigned int min_axis=0;
//     for(unsigned int i=0; i<3; i++)
//         if(un(i) > min_axis_size){
//             min_axis_size = un(i);
//             min_axis = i;
//         }
        
    
}



template<int dim> inline
void QXFEMFactory::refine_level(unsigned int n_simplices_to_refine)
{
    unsigned int new_level_offset = simplices_.size();
    unsigned int subsim_per_simplex = 1 << dim;
    simplices_.reserve(simplices_.capacity() + subsim_per_simplex * n_simplices_to_refine);
    
    for(unsigned int i = level_offset_; i<simplices_.size(); i++) {
        AuxSimplex &s = simplices_[i];
        if(s.refine) {
            refine_simplex<dim>(s);
            s.active = false;
        }
    }
    level_offset_ = new_level_offset;
}


template<int dim> inline
void QXFEMFactory::refine_simplex(const AuxSimplex& aux_simplex)
{
    static const unsigned int n_subelements = 1 << dim,  //2^dim
                              n_old_nodes = RefElement<dim>::n_nodes,
                              n_new_nodes = RefElement<dim>::n_lines; // new points are in the center of lines
    
    //TODO: make this RefElement dependent
    static const std::vector<unsigned int> conn[] = {
        {}, //0D
        
        {0, 2,
         2, 1},
         
        {0, 3, 4,
         3, 1, 5,
         4, 5, 2,
         3, 5, 4},
         
        {0, 4, 5, 6,
         4, 1, 7, 8,
         5, 7, 2, 9,
         6, 8, 9, 3,
         4, 6, 8, 9,
         4, 7, 8, 9,
         4, 6, 5, 9,
         4, 5, 7, 9}
    }; 
//     DebugOut().fmt("level = {}, {}\n", aux_simplex.level, max_refinement_level_);
 
    ASSERT_DBG(dim == aux_simplex.nodes.size()-1);
    
    // auxiliary vectors
    std::vector<Point> nodes = aux_simplex.nodes;
    nodes.reserve(n_old_nodes+n_new_nodes);

    // create new points on the simplex surface
    for(unsigned int e=0; e < n_new_nodes; e++)
    {
        Point p = nodes[RefElement<dim>::interact(Interaction<0,1>(e))[0]]
                  +nodes[RefElement<dim>::interact(Interaction<0,1>(e))[1]];
        nodes.push_back( p / 2.0);
    }
   
    
    for(unsigned int i=0; i < n_subelements; i++)
    {
        simplices_.push_back(AuxSimplex()); //create with active true, refine false
        AuxSimplex& s = simplices_.back();
        s.sing_id = aux_simplex.sing_id;
        
        s.nodes.resize(n_old_nodes);
        
        // over nodes
        for(unsigned int j=0; j < n_old_nodes; j++)
        {
            unsigned int conn_id = (n_old_nodes)*i + j;
            s.nodes[j] = nodes[conn[dim][conn_id]];
        }
    }
}


template<int dim, class S>
void QXFEMFactory::distribute_qpoints(std::vector<Point>& real_points,
                                      std::vector<double>& weights,
                                      const std::vector<S> & sing,
                                      double ele_measure)
{
    const unsigned int order = 3;
    QGauss<dim> gauss(order);
    const unsigned int n_gauss = gauss.size();
    std::vector<typename Space<dim+1>::Point> bary_points(n_gauss);
    
    // compute barycentric coordinates of gauss points
    for(unsigned int q=0; q<n_gauss; q++)
        bary_points[q] = RefElement<dim>::local_to_bary(gauss.point(q));
    
//     DebugOut() << "gauss size " << n_gauss << "\n";
    
    //TODO: compute number of active simplices during refinement
    unsigned int n_active_simplices = 0;
    for(unsigned int i=0; i<simplices_.size(); i++)
        if(simplices_[i].active) n_active_simplices++;
        
    real_points.reserve(n_gauss*n_active_simplices);
    weights.reserve(n_gauss*n_active_simplices);
    for(unsigned int i=0; i<simplices_.size(); i++)
    {
        const AuxSimplex& s = simplices_[i];
        if(s.active)
        {
            Point p;
            double measure = s.measure<dim>();
            bool check_qpoint_inside_sing = s.refine && (s.sing_id >=0) && (i>level_offset_);
            
            for(unsigned int q=0; q<n_gauss; q++)
            {
                p.zeros();
                //map gauss points to real
                for(unsigned int i=0; i<dim+1; i++)
                    p = p + bary_points[q][i]*s.nodes[i];
                
                //skip points inside singularity
                if(check_qpoint_inside_sing)   //simplex is intersecting singularity
                {
                    if(point_in_singularity<dim>(p, sing[s.sing_id])) continue;
                }
                
                real_points.push_back(p);
                
                // We want to have: sum JxW_q = |T|-|S_w|
                // which means for an affine mapping: sum w_q * |J| = |T| - |S_w|
                // therefore for T without singularity, it holds: sum w_q = |T| / |J|
                // which is for triangle equal 0.5, tetrahedron 0.1667
                double weight = gauss.weight(q) * measure / (ele_measure);
                weights.push_back(weight);
            }
        }
    }
    real_points.shrink_to_fit();
    weights.shrink_to_fit();
}

template<>
bool QXFEMFactory::point_in_singularity<2, QXFEMFactory::Singularity0DPtr>(const Point& p, Singularity0DPtr s){
    return s->geometry().point_in_ellipse(p);
}
template<>
bool QXFEMFactory::point_in_singularity<3, QXFEMFactory::Singularity1DPtr>(const Point& p, Singularity1DPtr s){
    return s->geometry().point_in_cylinder(p);
}

template<> double QXFEMFactory::AuxSimplex::measure<2>() const{
    ASSERT_DBG(nodes.size() == 3);
    Point s0 = nodes[1] - nodes[0],   // 0. edge of triangle
          s1 = nodes[2] - nodes[0];   // 1. edge of triangle
    return 0.5*arma::norm(arma::cross(s0,s1),2);
}

template<> double QXFEMFactory::AuxSimplex::measure<3>() const{
    ASSERT_DBG(nodes.size() == 4);
    Point s0 = nodes[1] - nodes[0],   // 0. edge
          s1 = nodes[2] - nodes[0],   // 1. edge
          s2 = nodes[3] - nodes[0];   // 2. edge
    return std::abs(arma::dot(arma::cross(s0,s1),s2)) / 6.0;
}
                
template<int dim>
void QXFEMFactory::map_real_to_unit_points(const std::vector<Point>& real_points,
                                           std::vector<typename Space<dim>::Point>& unit_points,
                                           ElementFullIter ele)
{
    ASSERT_DBG(ele->dim() == dim);
    unit_points.resize(real_points.size());
    
    MappingP1<dim,3> map;
    auto ele_map = map.element_map(*ele);
    
    for(unsigned int i=0; i<real_points.size(); i++)
    {
        unit_points[i] = RefElement<dim>::bary_to_local(map.project_real_to_unit(real_points[i],ele_map));
    }
}


template<int dim>
void QXFEMFactory::gnuplot_refinement(ElementFullIter ele,
                                      const string& output_dir,
                                      const QXFEM<dim,3>& quad,
                                      const std::vector<Singularity0DPtr> & sing)
{
 
    DBGCOUT("XFEM quadrature gnuplot print.\n");
    DBGVAR(output_dir);
    std::string fgnuplot_ref = "adaptive_integration_refinement_",
              fgnuplot_qpoints = "adaptive_integration_qpoints_",
              script_file = "g_script_adapt_",
              felements = "elements";
  
              fgnuplot_ref += std::to_string(ele->index()) + ".dat";
              fgnuplot_qpoints += std::to_string(ele->index()) + ".dat";
              script_file += std::to_string(ele->index()) + ".p";
        
    std::ofstream felements_file;
    felements_file.open (output_dir + felements, std::ios_base::app);
    if (felements_file.is_open()) 
    {
        for(unsigned int j=0; j<ele->n_nodes(); j++) {
            felements_file << ele->node[j]->getX() << " "
                << ele->node[j]->getY() << " " 
                << ele->node[j]->getZ() << "\n";
        }
        felements_file << ele->node[0]->getX() << " "
                << ele->node[0]->getY() << " " 
                << ele->node[0]->getZ() << "\n\n";
    }
    else 
    { 
        WarningOut() << "Coud not write refinement for gnuplot.\n";
    }
    felements_file.close();
    
    std::ofstream myfile1;
    myfile1.open (output_dir + fgnuplot_ref);
    if (myfile1.is_open()) {
        for (unsigned int i = 0; i < simplices_.size(); i++)
        {
            if(simplices_[i].active)
            for (unsigned int j = 0; j < dim+1; j++) 
            {
                Point &p = simplices_[i].nodes[j];
                myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n";
            }
            Point &p = simplices_[i].nodes[0];
            myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n\n";
        }
    }
    else 
    { 
        WarningOut() << "Coud not write refinement for gnuplot.\n";
    }
    myfile1.close();
    
    std::ofstream myfile2;
    myfile2.open (output_dir + fgnuplot_qpoints);
    if (myfile2.is_open()) 
    {
        for (unsigned int q = 0; q < quad.size(); q++) {
            const Point &p = quad.real_point(q);
            myfile2 << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        
    }
    else 
    { 
        WarningOut() << "Coud not write qpoints for gnuplot.\n";
    }
    myfile2.close();
    
    //creating command
    std::ostringstream strs;
    strs << "set terminal x11\n";
    strs << "set size ratio -1\n";
    strs << "set parametric\n";
    strs << "set xlabel 'X'\n";
    strs << "set ylabel 'Y'\n";
    strs << "set zlabel 'Z'\n";
       
    strs << "set style line 1 lt 2 lw 2 lc rgb 'blue'\n"
            << "set style line 2 lt 1 lw 2 lc rgb '#66A61E'\n";
    strs << "splot "; 
    
    for(unsigned int j = 0; j < sing.size(); j++)
    {
        ASSERT(dim==2);
        //print ellipse
        const CircleEllipseProjection& g = sing[j]->geometry();
        strs << "[t=0:2*pi] " 
            << g.center()[0] << " + " << g.ellipse_b()[0] << "*cos(t) + " <<  g.ellipse_a()[0] << "*sin(t), "
            << g.center()[1] << " + " << g.ellipse_b()[1] << "*cos(t) + " <<  g.ellipse_a()[1] << "*sin(t), "
            << g.center()[2] << " + " << g.ellipse_b()[2] << "*cos(t) + " <<  g.ellipse_a()[2] << "*sin(t) ls 1 ,\\\n";
        
        //print circle
        CircleEllipseProjection gg(g.center(),g.radius(),g.direction_vector(),g.direction_vector());
        strs << "[t=0:2*pi] " 
            << gg.center()[0] << " + " << gg.ellipse_b()[0] << "*cos(t) + " <<  gg.ellipse_a()[0] << "*sin(t), "
            << gg.center()[1] << " + " << gg.ellipse_b()[1] << "*cos(t) + " <<  gg.ellipse_a()[1] << "*sin(t), "
            << gg.center()[2] << " + " << gg.ellipse_b()[2] << "*cos(t) + " <<  gg.ellipse_a()[2] << "*sin(t) ls 1 ,\\\n";
    }
    
    strs << "'" << fgnuplot_qpoints << "' using 1:2:3 with points lc rgb 'light-blue' title 'quadrature points' ,\\\n";
    strs << "'" << fgnuplot_ref << "' using 1:2:3 with lines lc rgb 'red' title 'refinement'\n";
    
    //saving gnuplot script
    std::ofstream myfile3;
    myfile3.open (output_dir + script_file);
    if (myfile3.is_open()) 
    {
        // header
        myfile3 << "# Gnuplot script for printing adaptively refined element.\n" <<
                    "# Made by Pavel Exner.\n#\n" <<
                    "# Run the script in gnuplot:\n" <<
                    "# > load \"" << script_file << "\"\n#\n" <<
                    "# Data files used:\n" << 
                    "# " << fgnuplot_ref << "\n"
                    "# " << fgnuplot_qpoints << "\n" 
                    "#\n#" << std::endl;
        // script
        myfile3 << strs.str() << std::endl;
    }
    else 
    { 
        WarningOut() << "Coud not write gnuplot script.\n";
    }
    myfile3.close();
}

template<int dim>
void QXFEMFactory::gnuplot_refinement(ElementFullIter ele,
                                      const string& output_dir,
                                      const QXFEM<dim,3>& quad)
{
 
    DBGCOUT("XFEM quadrature gnuplot print.\n");
    DBGVAR(output_dir);
    std::string fgnuplot_ref = "adaptive_integration_refinement_",
              fgnuplot_qpoints = "adaptive_integration_qpoints_",
              script_file = "g_script_adapt_",
              felements = "elements";
  
              fgnuplot_ref += std::to_string(ele->index()) + ".dat";
              fgnuplot_qpoints += std::to_string(ele->index()) + ".dat";
              script_file += std::to_string(ele->index()) + ".p";
        
    std::ofstream felements_file;
    felements_file.open (output_dir + felements, std::ios_base::app);
    if (felements_file.is_open()) 
    {
        for(unsigned int j=0; j<ele->n_nodes(); j++) {
            felements_file << ele->node[j]->getX() << " "
                << ele->node[j]->getY() << " " 
                << ele->node[j]->getZ() << "\n";
        }
        felements_file << ele->node[0]->getX() << " "
                << ele->node[0]->getY() << " " 
                << ele->node[0]->getZ() << "\n\n";
    }
    else 
    { 
        WarningOut() << "Coud not write refinement for gnuplot.\n";
    }
    felements_file.close();
    
    std::ofstream myfile1;
    myfile1.open (output_dir + fgnuplot_ref);
    if (myfile1.is_open()) {
        for (unsigned int i = 0; i < simplices_.size(); i++)
        {
            if(simplices_[i].active)
//             if(simplices_[i].refine && (simplices_[i].sing_id >=0) && (i>level_offset_))
            {
                for (unsigned int j = 0; j < dim+1; j++) 
                {
                    Point &p = simplices_[i].nodes[j];
                    myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n";
                }
                Point &p = simplices_[i].nodes[0];
                myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n";
                p = simplices_[i].nodes[2];
                myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n";
                p = simplices_[i].nodes[1];
                myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n";
                p = simplices_[i].nodes[3];
                myfile1 << p[0] << " " << p[1] << " " << p[2] << "\n\n\n";
            }
        }
    }
    else 
    { 
        WarningOut() << "Coud not write refinement for gnuplot.\n";
    }
    myfile1.close();
    
    std::ofstream myfile2;
    myfile2.open (output_dir + fgnuplot_qpoints);
    if (myfile2.is_open()) 
    {
        for (unsigned int q = 0; q < quad.size(); q++) {
            const Point &p = quad.real_point(q);
            myfile2 << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        
    }
    else 
    { 
        WarningOut() << "Coud not write qpoints for gnuplot.\n";
    }
    myfile2.close();
    
    //creating command
    std::ostringstream strs;
    strs << "set terminal x11\n";
    strs << "set size ratio -1\n";
    strs << "set parametric\n";
    strs << "set xlabel 'X'\n";
    strs << "set ylabel 'Y'\n";
    strs << "set zlabel 'Z'\n";
       
    strs << "set style line 1 lt 2 lw 2 lc rgb 'blue'\n"
            << "set style line 2 lt 1 lw 2 lc rgb '#66A61E'\n";
    strs << "splot "; 
    
//     for(unsigned int j = 0; j < sing.size(); j++)
//     {
//         ASSERT(dim==2);
//         //print ellipse
//         const CircleEllipseProjection& g = sing[j]->geometry();
//         strs << "[t=0:2*pi] " 
//             << g.center()[0] << " + " << g.ellipse_b()[0] << "*cos(t) + " <<  g.ellipse_a()[0] << "*sin(t), "
//             << g.center()[1] << " + " << g.ellipse_b()[1] << "*cos(t) + " <<  g.ellipse_a()[1] << "*sin(t), "
//             << g.center()[2] << " + " << g.ellipse_b()[2] << "*cos(t) + " <<  g.ellipse_a()[2] << "*sin(t) ls 1 ,\\\n";
//         
//         //print circle
//         CircleEllipseProjection gg(g.center(),g.radius(),g.direction_vector(),g.direction_vector());
//         strs << "[t=0:2*pi] " 
//             << gg.center()[0] << " + " << gg.ellipse_b()[0] << "*cos(t) + " <<  gg.ellipse_a()[0] << "*sin(t), "
//             << gg.center()[1] << " + " << gg.ellipse_b()[1] << "*cos(t) + " <<  gg.ellipse_a()[1] << "*sin(t), "
//             << gg.center()[2] << " + " << gg.ellipse_b()[2] << "*cos(t) + " <<  gg.ellipse_a()[2] << "*sin(t) ls 1 ,\\\n";
//     }
//     
    strs << "'" << fgnuplot_qpoints << "' using 1:2:3 with points lc rgb 'light-blue' title 'quadrature points' ,\\\n";
    strs << "'" << fgnuplot_ref << "' using 1:2:3 with lines lc rgb 'red' title 'refinement',\\\n";
    strs << "'q_points.dat' using 1:2:3 with points lc rgb 'green' title 'q_points'\n";
    
    //saving gnuplot script
    std::ofstream myfile3;
    myfile3.open (output_dir + script_file);
    if (myfile3.is_open()) 
    {
        // header
        myfile3 << "# Gnuplot script for printing adaptively refined element.\n" <<
                    "# Made by Pavel Exner.\n#\n" <<
                    "# Run the script in gnuplot:\n" <<
                    "# > load \"" << script_file << "\"\n#\n" <<
                    "# Data files used:\n" << 
                    "# " << fgnuplot_ref << "\n"
                    "# " << fgnuplot_qpoints << "\n" 
                    "#\n#" << std::endl;
        // script
        myfile3 << strs.str() << std::endl;
    }
    else 
    { 
        WarningOut() << "Coud not write gnuplot script.\n";
    }
    myfile3.close();
}

// explicit instances:
template void QXFEMFactory::gnuplot_refinement<2>(ElementFullIter ele,
                                                  const string& output_dir,
                                                  const QXFEM<2,3>& quad,
                                                  const std::vector<Singularity0DPtr> & sing);
template void QXFEMFactory::gnuplot_refinement<3>(ElementFullIter ele,
                                                  const string& output_dir,
                                                  const QXFEM<3,3>& quad,
                                                  const std::vector<Singularity0DPtr> & sing);
template void QXFEMFactory::gnuplot_refinement<2>(ElementFullIter ele,
                                                  const string& output_dir,
                                                  const QXFEM<2,3>& quad);
template void QXFEMFactory::gnuplot_refinement<3>(ElementFullIter ele,
                                                  const string& output_dir,
                                                  const QXFEM<3,3>& quad);

template void QXFEMFactory::refine_level<2>(unsigned int n_simplices_to_refine);
template void QXFEMFactory::refine_level<3>(unsigned int n_simplices_to_refine);

template void QXFEMFactory::refine_simplex<2>(const AuxSimplex& aux_simplex);
template void QXFEMFactory::refine_simplex<3>(const AuxSimplex& aux_simplex);
