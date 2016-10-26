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

template<int dim, int spacedim>
const double QXFEMFactory<dim,spacedim>::distance_criteria_factor_ = 0.25;

template<int dim, int spacedim>
void QXFEMFactory< dim, spacedim >::clear()
{
    simplices_.clear();
    level_ = 0;
    level_offset_ = 0;
}

template<>
std::shared_ptr< QXFEM< 2, 2 > > QXFEMFactory<2,2>::create_singular(
                                                const std::vector<SingularityPtr> & sing,
                                                ElementFullIter ele)
{
    clear();
    std::shared_ptr<QXFEM<2,2>> qxfem = make_shared<QXFEM<2,2>>();
    ASSERT(0).error("Not implemented!");
    return qxfem;
}

template<>
std::shared_ptr< QXFEM< 2,3 > > QXFEMFactory<2,3>::create_singular(
                                                const std::vector<SingularityPtr> & sing,
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
    
    // WORKS only for single well
    // project the simplex into the plane of the circle
    sing[0]->geometry().project_to_circle_plane(s.nodes);
    
    s.active = true;
    simplices_.push_back(s);
    
    for(unsigned int i=0; i < max_level_; i++)
    {
        DebugOut() << "QXFEM level: " << i << "\n";
        unsigned int n_to_refine = refine_edge(sing);
        if(n_to_refine == 0) break; // end refinement
        
        refine_level(n_to_refine);
    }
    refine_edge(sing);
    
    // project the simplices back to the plane of the ellipse
    for(AuxSimplex& simp : simplices_){
        sing[0]->geometry().project_to_ellipse_plane(simp.nodes);
    }
    
    distribute_qpoints(qxfem->real_points_, qxfem->weights, sing, ele->measure());
    map_real_to_unit_points(qxfem->real_points_, qxfem->quadrature_points, ele);
    DebugOut().fmt("n_real_qpoints {} {}\n",qxfem->real_points_.size(), qxfem->quadrature_points.size());
    
    return qxfem;
}


template<int dim, int spacedim>
void QXFEMFactory<dim,spacedim>::refine_level(unsigned int n_simplices_to_refine)
{
    unsigned int new_level_offset = simplices_.size();
    unsigned int subsim_per_simplex = 1 << dim;
    simplices_.reserve(simplices_.capacity() + subsim_per_simplex * n_simplices_to_refine);
    
    for(unsigned int i = level_offset_; i<simplices_.size(); i++) {
        AuxSimplex &s = simplices_[i];
        if(s.refine) {
            refine_simplex(s);
            s.active = false;
        }
    }
    level_offset_ = new_level_offset;
}


template<int dim, int spacedim>
unsigned int QXFEMFactory<dim,spacedim>::refine_edge(const std::vector<SingularityPtr> & sing)
{
    ASSERT_DBG(dim == 2);
    ASSERT_DBG( (spacedim == 2) || (spacedim == 3));
    
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
            double max_h;
            int res = simplex_sigularity_intersection(*sing[j],s, distance_sqr, max_h);
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
                if( max_h > distance_criteria_factor_ * rmin)
                //if( max_h > distance_criteria_factor_ * distance_sqr)
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

template<int dim, int spacedim>
int QXFEMFactory<dim,spacedim>::simplex_sigularity_intersection(const Singularity0D<spacedim>& w,
                                                                const AuxSimplex& s,
                                                                double& distance_sqr,
                                                                double& max_h)
{
    ASSERT_DBG(dim == 2);
    ASSERT_DBG( (spacedim == 2) || (spacedim == 3));
    
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
            
            max_h = std::max(std::max(d00,d11),d22);
            
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


template<int dim, int spacedim>
void QXFEMFactory<dim,spacedim>::refine_simplex(const AuxSimplex& aux_simplex)
{
    static const unsigned int n_subelements = 1 << dim,  //2^dim
                              n_old_nodes = RefElement<dim>::n_nodes,
                              n_new_nodes = RefElement<dim>::n_lines; // new points are in the center of lines
                       
    static const std::vector<unsigned int> conn[] = {
        {}, //0D
        
        {0, 2,
         2, 1},
         
        {0, 3, 4,
         3, 1, 5,
         4, 5, 2,
         3, 5, 4},
         
        {0, 4, 5, 7,
         4, 1, 6, 8,
         5, 6, 2, 9,
         7, 8, 9, 3,
         4, 7, 8, 9,
         4, 6, 8, 9,
         4, 7, 5, 9,
         4, 5, 6, 9}
    }; 
//     DebugOut().fmt("level = {}, {}\n", aux_simplex.level, max_refinement_level_);
 
    ASSERT_DBG(dim == aux_simplex.nodes.size()-1);
    
    
    
    // auxiliary vectors
    std::vector<Point> nodes = aux_simplex.nodes;
    nodes.reserve(n_old_nodes+n_new_nodes);

    // create new points on the simplex surface
    for(unsigned int e=0; e < n_new_nodes; e++)
    {
        Point p = nodes[RefElement<dim>::template interact<0,1>(e)[0]]
                  +nodes[RefElement<dim>::template interact<0,1>(e)[1]];
        nodes.push_back( p / 2.0);
        //nodes.back().print();
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

template<int dim, int spacedim>
void QXFEMFactory< dim, spacedim >::distribute_qpoints(std::vector<Point>& real_points,
                                                       std::vector<double>& weights,
                                                       const std::vector<SingularityPtr> & sing,
                                                       double ele_measure)
{
    const unsigned int order = 3;
    QGauss<dim> gauss(order);
    const unsigned int n_gauss = gauss.size();
    std::vector<typename Space<dim+1>::Point> bary_points(n_gauss);
    
    // compute barycentric coordinates of gauss points
    for(unsigned int q=0; q<n_gauss; q++)
    {
        // complementary barycentric coord
        double complement = 1.0;
        for(unsigned int i=0; i<dim; i++) {
            complement = complement - gauss.point(q)[i];
            bary_points[q][i] = gauss.point(q)[i];
        }
        bary_points[q][dim] = complement;
    }
//     DebugOut() << "gauss size " << n_gauss << "\n";
    
    //TODO: compute number of active simplices during refinement
    unsigned int n_active_simplices = 0;
    for(unsigned int i=0; i<simplices_.size(); i++)
        if(simplices_[i].active) n_active_simplices++;
        
    real_points.reserve(order*n_active_simplices);
    weights.reserve(order*n_active_simplices);
    for(unsigned int i=0; i<simplices_.size(); i++)
    {
        const AuxSimplex& s = simplices_[i];
        if(s.active)
        {
            Point p;
            Point s0 = s.nodes[1] - s.nodes[0],   // 0. edge of triangle
                  s1 = s.nodes[2] - s.nodes[0];   // 1. edge of triangle
            double area = arma::norm(arma::cross(s0,s1),2);
            
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
                    if(sing[s.sing_id]->geometry().point_in_ellipse(p)) continue;
                }
                
                real_points.push_back(p);
                
                // We want to have: sum JxW_q = |T|-|S_w|
                // which means for an affine mapping: sum w_q * |J| = |T| - |S_w|
                // therefore for T without singularity, it holds: sum w_q = |T| / |J|
                // which is for 2D equal 0.5
                double weight = gauss.weight(q) * area / (2*ele_measure);
                weights.push_back(weight);
            }
        }
    }
    real_points.shrink_to_fit();
    weights.shrink_to_fit();
}

template<int dim, int spacedim>
void QXFEMFactory<dim,spacedim>::map_real_to_unit_points(const std::vector<Point>& real_points,
                                                         std::vector< typename Space< dim >::Point >& unit_points,
                                                         ElementFullIter ele)
{
    ASSERT(dim == 2);
    ASSERT_DBG( (spacedim == 2) || (spacedim == 3));
    ASSERT_DBG(ele->dim() == dim);
    
    unit_points.resize(real_points.size());
    
    // using barycentric coordinates
    // http://www.blackpawn.com/texts/pointinpoly/
    Point v0 = ele->node[1]->point() - ele->node[0]->point(),   // 0. edge of triangle
          v1 = ele->node[2]->point() - ele->node[0]->point();   // 1. edge of triangle
    // v2 = P - A
    // Point in triangle is defined: 
    // v2 = s*v0 + t*v1;
    // make 2 equations for s and t:
    // v2.v0 = s*v0.v0 + t*v1.v0;
    // v2.v1 = + s*v0.v1 + t*v1.v1;
    
    // dot products:
    double
        d00 = arma::dot(v0,v0),
        d01 = arma::dot(v0,v1),
        d11 = arma::dot(v1,v1);
    double invdenom = 1.0/ (d00*d11 - d01*d01);
    
    for(unsigned int i=0; i<real_points.size(); i++)
    {
        const Point& p = real_points[i];
        //map real to reference element
        Point v2 = p - simplices_[0].nodes[0];
        double d02 = arma::dot(v0,v2), 
                d12 = arma::dot(v1,v2);
        double t0 = (d02*d11 - d12*d01) * invdenom,
                t1 = (d12*d00 - d02*d01) * invdenom;
            
        unit_points[i] = {t0, t1};
    }
}


template<int dim, int spacedim>
void QXFEMFactory<dim,spacedim>::gnuplot_refinement(ElementFullIter ele,
                                                    const string& output_dir,
                                                    const QXFEM<dim,spacedim>& quad,
                                                    const std::vector<SingularityPtr> & sing)
{
    
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


// explicit instances:
template class QXFEMFactory<2,2>;
template class QXFEMFactory<2,3>;
