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
 * @file    tracing_algorithm.cc
 * @ingroup intersection
 * @brief   Implementation of tracing algorithms for polygonal intersections.
 */

#include "tracing_algorithm.hh"
#include "intersection_point_aux.hh"
#include "intersection_aux.hh"

#include "mesh/ref_element.hh"
#include <algorithm>

#include "system/sys_profiler.hh"

namespace computeintersection{
/*
void Tracing::trace_polygon(IntersectionAux<2,3> &p){
    
//     DebugOut().fmt("{} intersections:\n",p.size());
//     for(IntersectionPointAux<2,3> &ip : p.points())
//         DebugOut() << ip;
    
    // avoid tracing (none is needed) if the intersection is just single point
    if(p.points().size() < 2) {
        DebugOut() << "Less than 2 IPs, no tracing.\n";
        return;
    }
    
    if(p.is_pathologic()) trace_polygon_convex_hull(p);
    else trace_polygon_opt(p);
};

void Tracing::trace_polygon_opt(IntersectionAux<2,3> &p){

    START_TIMER("CI trace opt");
    ASSERT_DBG(!p.is_pathologic());

    const unsigned int tt_size = 7;
    // Tracing table - 7 (4 sides of tetrahedron + 3 edges of triangle) X 2
    // trace_table_x[i][j]: 
    // i - index of next row where we continue tracing
    // j - index of IP in the original vector of IPs
    vector<vector<int>> trace_table_x(tt_size);
    for(int i = 0; i < (int)tt_size;i++){
        trace_table_x[i].resize(2,-1);
    }
    
    // index of the first set row
    unsigned int first_row_index = tt_size;
    
    // go through all intersection points (vertices of polygon)
    for(unsigned int i = 0; i < p.points().size();i++){
        IntersectionPointAux<2,3> ip = p.points()[i];
        
        unsigned int row, object_index;
        
        switch(ip.dim_A())
        {
            case 0: // if the ip is at the vertex of triangle -> intersection E-E
            {
                // IP is the vertex of triangle,
                // the directions of triangle edges determines the directions of polygon edges
                unsigned int vertex_index = ip.idx_A();
                DebugOut() << "Intersection type E-E. (Vertex_index: " << vertex_index << ")\n";
                
                row = 4 + RefElement<2>::interact(Interaction<1,0>(vertex_index))[0];
                object_index = 4 + RefElement<2>::interact(Interaction<1,0>(vertex_index))[1];
                DebugOut().fmt("E-E: row: {}, to edge: {}, ip: {} \n",row, object_index, i);

            } break;
            case 1: // if the ip is on the side of triangle -> intersection S-E or E-S
            {

                DebugOut().fmt("Intersection type E{} - S{}, ori {}.\n",ip.idx_A(),ip.idx_B(),int(ip.orientation()));
                // IP is not a vertex of triangle -> cannot be E-E
                /* directions of polygon edges depends on three aspects:
                 * - normal orientation of the tetrahedron side (OUT...0, IN...1)
                 * - orientation of the edge of the triangle (according to direction of other edges of triangle or opposite)
                 * - orientation of the intersection point ... Plucker product (P>0...1, P<0...0)
                 * 
                 * we can obtain this XOR table: XOR(XOR(S,P),E)
                 *  S   P   E   RES
                 *  0   0   0   0
                 *  0   0   1   1
                 *  0   1   0   1
                 *  0   1   1   0
                 *  1   0   0   1
                 *  1   0   1   0
                 *  1   1   0   0
                 *  1   1   1   1
                 * 
                 * -> this can be expressed as sum(S,P,E)%2
                 * 
                 * if RES == 1 -> intersection direction is E-S
                 * if RES == 0 -> intersection direction is S-E (tetrahedron Side -> triangle Edge)
                 *//*
                ASSERT_DBG( ip.dim_B() == 2 );
                unsigned int sign = (unsigned int)(ip.orientation());
                ASSERT_DBG( sign < 2);
                unsigned int j = (RefElement<3>::normal_orientation(ip.idx_B()) +
                                  RefElement<2>::normal_orientation(ip.idx_A()) +
                                  sign
                                 ) % 2;
                if(j == 0){
                    // direction: Side -> IP -> Edge
                    row = ip.idx_B();
                    object_index = 4 + ip.idx_A();

                    DebugOut().fmt("S-E: row: {}, to edge: {}, ip: {} \n",row, object_index, i);
                }else{
                    // direction: Edge -> IP -> Side
                    row = 4 + ip.idx_A();
                    object_index = ip.idx_B();

                    DebugOut().fmt("E-S: row: {}, to side: {}, ip: {} \n",row, object_index, i);
                }
            } break;
            case 2: // type S-S, IP is on the edge between two sides
            {
                ASSERT_DBG( ip.dim_B() == 1 );
                unsigned int tetrahedron_line = ip.idx_B();
                unsigned int first_side = (unsigned int)(ip.orientation());
                ASSERT_DBG( first_side < 2);

                DebugOut().fmt("Intersection type S-S, ori {}.\n", first_side);
                /** here the idx_B contains number of edge of the tetrahedron where the IP lies
                 * sides: let edge be oriented up and let us see the tetrahedron from outside ->
                 *  then the right side [0] is in and left side [1] is out
                 *  - this is determined by IP orientation
                 *//*

                DebugOut().VarFmt(ip.idx_B()).VarFmt(first_side);
                row = RefElement<3>::interact(Interaction<2,1>(tetrahedron_line))[1-first_side];
                object_index = RefElement<3>::interact(Interaction<2,1>(tetrahedron_line))[first_side];

                DebugOut().fmt("S-S: row: {}, to side: {}, ip: {} \n",row, object_index, i);
            } break;
            default:
                ASSERT_DBG(0).error("Unsupported dimension of intersection object A.");
                row = object_index = 0; // suppress compilation warnings
              break;
        }
        
        // determine the first row set
        if(row < first_row_index) first_row_index = row;
        
        trace_table_x[row][0] = object_index;
        trace_table_x[row][1] = i;
    }
    
//     DebugOut() << "Trace table:\n";
//     for(unsigned int i = 0; i < 7;i++){
//         DebugOut() << "i=" << i << "  ";
//         for(unsigned int j = 0; j < 2;j++)
//             DebugOut() << trace_table_x[i][j] << "  ";
//         DebugOut() << endl;
//     }
    ASSERT_LT_DBG(first_row_index, tt_size);
    
    
    // traced IPs (reordered)
    std::vector<IntersectionPointAux<2,3>> new_points;
    new_points.reserve(p.points().size());
    
    // start polygon with the IP at the first row set
    new_points.push_back(p.points()[trace_table_x[first_row_index][1]]);
    // determine next row in the trace table
    unsigned int next_row = trace_table_x[first_row_index][0];

    // jump from row to row until we get back to the first row 
    // (i.e. go through polygonal vertices until we get back to starting point)
    while(next_row != first_row_index){
        DebugOut().VarFmt(next_row);
        unsigned int i_ip_orig = (unsigned int)trace_table_x[next_row][1];
        ASSERT_LT_DBG(i_ip_orig, p.points().size());
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // one ugly awful HACK:
        // if the optimal tracing fails, then call convex hull tracing
        if(i_ip_orig >= p.points().size())
        {
            trace_polygon_convex_hull(p);
            break;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        new_points.push_back(p.points()[i_ip_orig]);
        next_row = trace_table_x[next_row][0];
    }

    p.points() = new_points;

    END_TIMER("CI trace opt");
};

/*
void Tracing::trace_polygon_convex_hull(IntersectionAux<2,3> &p){

    START_TIMER("CI trace convex hull");
    DebugOut() << "convex hull tracing\n";
    
    // sort IPs using IP's operator <
    std::sort(p.points().begin(),p.points().end());

    //TODO: think about removing?? at this point we have already found IPs, no need to remove them...
    // Remove duplicit IPs, use geometry epsilon for comparison
    for(unsigned int i = 0; i < p.points().size()-1; i++){
        if((fabs(p.points()[i].local_bcoords_A()[0] - p.points()[i+1].local_bcoords_A()[0]) < geometry_epsilon) &&
           (fabs(p.points()[i].local_bcoords_A()[1] - p.points()[i+1].local_bcoords_A()[1]) < geometry_epsilon)){
            p.points().erase(p.points().begin()+i);
        }
    }

    // check also the first and last IP
    if(p.points().size() > 1 && fabs(p.points()[0].local_bcoords_A()[0] - p.points()[p.points().size()-1].local_bcoords_A()[0]) < geometry_epsilon &&
            fabs(p.points()[0].local_bcoords_A()[1] - p.points()[p.points().size()-1].local_bcoords_A()[1]) < geometry_epsilon){
        p.points().erase(p.points().end());
    }

    DebugOut() << "Convex hull - number of points: " << p.points().size() << "\n";
    
    int n = p.points().size(), k = 0;
    if(n > 2){
    
        std::vector<IntersectionPointAux<2,3>> H(2*n);

        for(int i = 0; i < n; ++i){
            while(k >= 2 && convex_hull_cross(H[k-2], H[k-1], p.points()[i]) <= 0) k--;
            H[k++] = p.points()[i];
        }

        for(int i = n-2, t = k+1;i>=0;i--){
            while(k >= t && (convex_hull_cross(H[k-2], H[k-1], p.points()[i])) <= 0) k--;
            H[k++] = p.points()[i];
        }

        H.resize(k-1);
        p.points() = H;
    }

    END_TIMER("CI trace convex hull");
};

double Tracing::convex_hull_cross(const IntersectionPointAux<2,3> &O,
        const IntersectionPointAux<2,3> &A,
        const IntersectionPointAux<2,3> &B){
    return ((A.local_bcoords_A()[1]-O.local_bcoords_A()[1])*(B.local_bcoords_A()[2]-O.local_bcoords_A()[2])
            -(A.local_bcoords_A()[2]-O.local_bcoords_A()[2])*(B.local_bcoords_A()[1]-O.local_bcoords_A()[1]));
}

int Tracing::side_content_prolongation(IntersectionAux<2,3> &p){

    vector<unsigned int> side_field(4,0);

    // compute number of IPs lying on each tetrahedron face
    for(unsigned int i = 0; i < p.points().size();i++){

        for(unsigned int j = 0; j < 4;j++){
            if(fabs((double)p.points()[i].local_bcoords_B()[j]) < geometry_epsilon){
                side_field[j]++;
            }
        }
    }

    
    for(unsigned int i = 0;i< 4;i++){
        if(side_field[i] > 2){
            return i;
        }
    }
    return -1;

}*/

} // END NAMESPACE
