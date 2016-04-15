/*
 * IntersectionPolygon.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "trace_algorithm.h"
#include "intersectionpoint.h"
#include "intersectionaux.h"

#include "mesh/ref_element.hh"
#include <algorithm>

namespace computeintersection{
    
void Tracing::trace_polygon(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &p){
    
    if(p.is_pathologic()) trace_polygon_convex_hull(prolongation_table, p);
    else trace_polygon_opt(prolongation_table, p);
};

void Tracing::trace_polygon_opt(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &p){

    ASSERT(!p.is_pathologic(), "Cannot call this polygonal tracing in pathologic case.");
    
    prolongation_table.clear();
    
    //TODO: this is checked outside in inspect_elements; assert??
    // avoid tracing (none is needed) if the intersection is just single point
    if(p.points().size() < 2) return;

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

//     for(unsigned int i = 0; i < i_points_.size();i++)
//         i_points_[i].print();    
    
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
                DBGMSG("Intersection type E-E. (Vertex_index: %d)\n", vertex_index);
                
                row = 4 + RefElement<2>::interact<1,0>(vertex_index)[0];
                object_index = 4 + RefElement<2>::interact<1,0>(vertex_index)[1];
                DBGMSG("E-E: row: %d, to edge: %d, ip: %d \n",row, object_index, i);

            } break;
            case 1: // if the ip is on the side of triangle -> intersection S-E or E-S
            {
                DBGMSG("Intersection type E%d - S%d, ori %d.\n",ip.idx_A(),ip.idx_B(),ip.orientation());
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
                 * if RES == 1 -> intersection direction is S-E (tetrahedron Side -> triangle Edge)
                 * if RES == 0 -> intersection direction is E-S
                 */
                unsigned int j = (RefElement<3>::normal_orientation(ip.idx_B()) + 
                                  RefElement<2>::normal_orientation(ip.idx_A()) +
                                  ip.orientation()
                                 ) % 2;
                if(j == 1){
                    // direction: Side -> IP -> Edge
                    row = ip.idx_B();
                    object_index = 4 + ip.idx_A();

                    DBGMSG("S-E: row: %d, to edge: %d, ip: %d \n",row, object_index, i);
                }else{
                    // direction: Edge -> IP -> Side
                    row = 4 + ip.idx_A();
                    object_index = ip.idx_B();

                    DBGMSG("E-S: row: %d, to side: %d, ip: %d \n",row, object_index, i);
                }
            } break;
            case 2: // type S-S, IP is on the edge between two sides
            {
                DBGMSG("Intersection type S-S, ori %d.\n", ip.orientation());
                /** here the idx_B contains number of edge of the tetrahedron where the IP lies
                 * sides: let edge be oriented up and let us see the tetrahedron from outside ->
                 *  then the right side [0] is in and left side [1] is out
                 *  - this is determined by IP orientation
                 */
                unsigned int tetrahedron_line = ip.idx_B();
                row = RefElement<3>::interact<2,1>(tetrahedron_line)[1-ip.orientation()];
                object_index = RefElement<3>::interact<2,1>(tetrahedron_line)[ip.orientation()];

                DBGMSG("S-S: row: %d, to side: %d, ip: %d \n",row, object_index, i);
            } break;
            default:
                ASSERT(0,"Unsupported dimension of intersection object A.");
                row = object_index = 0; // suppress compilation warnings
              break;
        }
        
        // determine the first row set
        if(row < first_row_index) first_row_index = row;
        
        trace_table_x[row][0] = object_index;
        trace_table_x[row][1] = i;
    }
    
//     DBGMSG("Trace table:\n");
//     for(unsigned int i = 0; i < 7;i++){
//         cout << "i=" << i << "  ";
//         for(unsigned int j = 0; j < 2;j++)
//             std::cout << trace_table_x[i][j] << "  ";
//         cout << endl;
//     }
    ASSERT_LESS(first_row_index,tt_size);
    
    
    // traced IPs (reordered)
    std::vector<IntersectionPointAux<2,3>> new_points;
    new_points.reserve(p.points().size());
    
    // start polygon with the IP at the first row set
    new_points.push_back(p.points()[trace_table_x[first_row_index][1]]);
    // determine next row in the trace table
    unsigned int next_row = trace_table_x[first_row_index][0];
    
    // Fill also the prolongation table
    // It contains indices of tetrahedron sides and triangle edges (indices in the tracing table)
    prolongation_table.reserve(p.points().size());
    prolongation_table.push_back((unsigned int)trace_table_x[first_row_index][0]);

    // jump from row to row until we get back to the first row 
    // (i.e. go through polygonal vertices until we get back to starting point)
    while(next_row != first_row_index){
        DBGMSG("next_row = %d\n",next_row);
        unsigned int i_ip_orig = (unsigned int)trace_table_x[next_row][1];
        ASSERT_LESS(i_ip_orig, p.points().size());
        new_points.push_back(p.points()[i_ip_orig]);
        prolongation_table.push_back((unsigned int)trace_table_x[next_row][0]);
        next_row = trace_table_x[next_row][0];
    }

    p.points() = new_points;

};

void Tracing::trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table, IntersectionAux<2,3> &p){

    DBGMSG("convex hull tracing\n");
    //TODO: this is checked outside in inspect_elements; assert??
    // skip tracing if not enough IPs
    if(p.points().size() <= 1) return;

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


    int n = p.points().size(), k = 0;
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


    // Filling prolongation table
    if(p.points().size() > 1){

        // Prolongation - finding same zero bary coord/coords except
        // side with content intersect
        unsigned int size = p.points().size() == 2 ? 1 : p.points().size();
        int forbidden_side = side_content_prolongation(p);

        for(unsigned int i = 0; i < size;i++){
            unsigned int ii = (i+1)%p.points().size();

            for(unsigned int j = 0; j < 3;j++){
                if(fabs((double)p.points()[i].local_bcoords_A()[j]) < geometry_epsilon && fabs((double)p.points()[ii].local_bcoords_A()[j]) < geometry_epsilon){
                    prolongation_table.push_back(6-j);
                }
            }
            for(unsigned int k = 0; k < 4;k++){
                if((int)k != forbidden_side && fabs((double)p.points()[i].local_bcoords_B()[k]) < geometry_epsilon && fabs((double)p.points()[ii].local_bcoords_B()[k]) < geometry_epsilon){
                    prolongation_table.push_back(3-k);
                }
            }
        }
    }

};

double Tracing::convex_hull_cross(const IntersectionPointAux<2,3> &O,
        const IntersectionPointAux<2,3> &A,
        const IntersectionPointAux<2,3> &B){
    return ((A.local_bcoords_A()[1]-O.local_bcoords_A()[1])*(B.local_bcoords_A()[2]-O.local_bcoords_A()[2])
            -(A.local_bcoords_A()[2]-O.local_bcoords_A()[2])*(B.local_bcoords_A()[1]-O.local_bcoords_A()[1]));
}

int Tracing::side_content_prolongation(IntersectionAux<2,3> &p){

    vector<unsigned int> side_field(4,0);

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

}

} // END NAMESPACE
