/*
 * IntersectionPolygon.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "intersectionpolygon.h"
#include "intersectionpoint.h"
#include "mesh/ref_element.hh"
#include <deal.II/base/geometry_info.h>
#include <algorithm>

namespace computeintersection{
IntersectionPolygon::IntersectionPolygon() {};

IntersectionPolygon::~IntersectionPolygon() {
	// TODO Auto-generated destructor stub
}

IntersectionPolygon::IntersectionPolygon(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){

	pathologic_ = false;
	i_points_.reserve(4);
};

void IntersectionPolygon::add_ipoint(const IntersectionPoint<2,3> &InPoint){
	if(InPoint.is_pathologic()){
		pathologic_ = true;
	}
	i_points_.push_back(InPoint);
};

void IntersectionPolygon::trace_polygon(std::vector<unsigned int> &prolongation_table){
    
    if(pathologic_) trace_polygon_convex_hull(prolongation_table);
    else trace_polygon_opt(prolongation_table);
};

void IntersectionPolygon::trace_polygon_opt(std::vector<unsigned int> &prolongation_table){

    ASSERT(!pathologic_, "Cannot call this polygonal tracing in pathologic case.");
    
    //TODO: this is checked outside in inspect_elements; assert??
    // avoid tracing (none is needed) if the intersection is just single point
	if(i_points_.size() < 2) return;

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
	for(unsigned int i = 0; i < i_points_.size();i++){
        IntersectionPoint<2,3> ip = i_points_[i];
        
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
    std::vector<IntersectionPoint<2,3>> new_points;
    new_points.reserve(i_points_.size());
    
    // start polygon with the IP at the first row set
    new_points.push_back(i_points_[trace_table_x[first_row_index][1]]);
    // determine next row in the trace table
    unsigned int next_row = trace_table_x[first_row_index][0];
    
    // Fill also the prolongation table
    // It contains indices of tetrahedron sides and triangle edges (indices in the tracing table)
    prolongation_table.reserve(i_points_.size());
    prolongation_table.push_back((unsigned int)trace_table_x[first_row_index][0]);

    // jump from row to row until we get back to the first row 
    // (i.e. go through polygonal vertices until we get back to starting point)
    while(next_row != first_row_index){
        DBGMSG("next_row = %d\n",next_row);
        unsigned int i_ip_orig = (unsigned int)trace_table_x[next_row][1];
        ASSERT_LESS(i_ip_orig, i_points_.size());
        new_points.push_back(i_points_[i_ip_orig]);
        prolongation_table.push_back((unsigned int)trace_table_x[next_row][0]);
        next_row = trace_table_x[next_row][0];
    }

	i_points_ = new_points;

};

void IntersectionPolygon::trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table){

    //TODO: this is checked outside in inspect_elements; assert??
    // skip tracing if not enough IPs
    if(i_points_.size() <= 1) return;

    // sort IPs using IP's operator <
	std::sort(i_points_.begin(),i_points_.end());

    //TODO: think about removing?? at this point we have already found IPs, no need to remove them...
	// Remove duplicit IPs, use geometry epsilon for comparison
	for(unsigned int i = 0; i < i_points_.size()-1; i++){
		if((fabs(i_points_[i].local_bcoords_A()[0] - i_points_[i+1].local_bcoords_A()[0]) < geometry_epsilon) &&
		   (fabs(i_points_[i].local_bcoords_A()[1] - i_points_[i+1].local_bcoords_A()[1]) < geometry_epsilon)){
			i_points_.erase(i_points_.begin()+i);
		}
	}

	// check also the first and last IP
	if(i_points_.size() > 1 && fabs(i_points_[0].local_bcoords_A()[0] - i_points_[i_points_.size()-1].local_bcoords_A()[0]) < geometry_epsilon &&
			fabs(i_points_[0].local_bcoords_A()[1] - i_points_[i_points_.size()-1].local_bcoords_A()[1]) < geometry_epsilon){
		i_points_.erase(i_points_.end());
	}


	int n = i_points_.size(), k = 0;
	std::vector<IntersectionPoint<2,3>> H(2*n);

	for(int i = 0; i < n; ++i){
		while(k >= 2 && convex_hull_cross(H[k-2], H[k-1], i_points_[i]) <= 0) k--;
		H[k++] = i_points_[i];
	}

	for(int i = n-2, t = k+1;i>=0;i--){
		while(k >= t && (convex_hull_cross(H[k-2], H[k-1], i_points_[i])) <= 0) k--;
		H[k++] = i_points_[i];
	}

	H.resize(k-1);
	i_points_ = H;


	// Filling prolongation table
	if(i_points_.size() > 1){

		// Prolongation - finding same zero bary coord/coords except
		// side with content intersect
		unsigned int size = i_points_.size() == 2 ? 1 : i_points_.size();
		int forbidden_side = side_content_prolongation();

		for(unsigned int i = 0; i < size;i++){
			unsigned int ii = (i+1)%i_points_.size();

			for(unsigned int j = 0; j < 3;j++){
				if(fabs((double)i_points_[i].local_bcoords_A()[j]) < geometry_epsilon && fabs((double)i_points_[ii].local_bcoords_A()[j]) < geometry_epsilon){
					prolongation_table.push_back(6-j);
				}
			}
			for(unsigned int k = 0; k < 4;k++){
				if((int)k != forbidden_side && fabs((double)i_points_[i].local_bcoords_B()[k]) < geometry_epsilon && fabs((double)i_points_[ii].local_bcoords_B()[k]) < geometry_epsilon){
					prolongation_table.push_back(3-k);
				}
			}
		}
	}

};

double IntersectionPolygon::convex_hull_cross(const IntersectionPoint<2,3> &O,
		const IntersectionPoint<2,3> &A,
		const IntersectionPoint<2,3> &B) const{
	return ((A.local_bcoords_A()[1]-O.local_bcoords_A()[1])*(B.local_bcoords_A()[2]-O.local_bcoords_A()[2])
			-(A.local_bcoords_A()[2]-O.local_bcoords_A()[2])*(B.local_bcoords_A()[1]-O.local_bcoords_A()[1]));
}

// int IntersectionPolygon::convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
//     		const IntersectionPoint<2,3> &B) const{
// 
// 	// Finding same zero 2D local coord
// 
// 	if(fabs((double)A.local_bcoords_A()[0]) < geometry_epsilon && fabs((double)B.local_bcoords_A()[0]) < epsilon){
// 		return 6;
// 	}else if(fabs((double)A.local_bcoords_A()[1]) < epsilon && fabs((double)B.local_bcoords_A()[1]) < epsilon){
// 		return 5;
// 	}else if(fabs((double)A.local_bcoords_A()[2]) < epsilon && fabs((double)B.local_bcoords_A()[2]) < epsilon){
// 		return 4;
// 	}else if(fabs((double)A.local_bcoords_B()[0]) < epsilon && fabs((double)B.local_bcoords_B()[0]) < epsilon){
// 		return 3;
// 	}else if(fabs((double)A.local_bcoords_B()[1]) < epsilon && fabs((double)B.local_bcoords_B()[1]) < epsilon){
// 		return 2;
// 	}else if(fabs((double)A.local_bcoords_B()[2]) < epsilon && fabs((double)B.local_bcoords_B()[2]) < epsilon){
// 		return 1;
// 	}else if(fabs((double)A.local_bcoords_B()[3]) < epsilon && fabs((double)B.local_bcoords_B()[3]) < epsilon){
// 		return 0;
// 	}
// 	// Todo: Nahradit Assertem
// 	/*cout << A.local_bcoords_A()[0] << ","
// 			<< A.local_bcoords_A()[1] << ","
// 			<< A.local_bcoords_A()[2] << ","
// 			<< A.local_bcoords_B()[0] << ","
// 			<< A.local_bcoords_B()[1] << ","
// 			<< A.local_bcoords_B()[2] << ","
// 			<< A.local_bcoords_B()[3] << "," << endl;
// 
// 
// 	cout << B.local_bcoords_A()[0] << ","
// 			<< B.local_bcoords_A()[1] << ","
// 			<< B.local_bcoords_A()[2] << ","
// 			<< B.local_bcoords_B()[0] << ","
// 			<< B.local_bcoords_B()[1] << ","
// 			<< B.local_bcoords_B()[2] << ","
// 			<< B.local_bcoords_B()[3] << ",";
// 	xprintf(Msg, "Mozny problem\n");*/
// 	return -1;
// }

int IntersectionPolygon::side_content_prolongation() const{

	vector<unsigned int> side_field(4,0);

	for(unsigned int i = 0; i < i_points_.size();i++){

		for(unsigned int j = 0; j < 4;j++){
			if(fabs((double)i_points_[i].local_bcoords_B()[j]) < geometry_epsilon){
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


double IntersectionPolygon::get_area() const{
	double subtotal = 0.0;
	for(unsigned int j = 2; j < i_points_.size();j++){
		//xprintf(Msg, "volani %d %d\n",j, i_points_.size());
		subtotal += fabs(i_points_[0].local_bcoords_A()(1)*(i_points_[j-1].local_bcoords_A()(2) - i_points_[j].local_bcoords_A()(2)) +
				 i_points_[j-1].local_bcoords_A()(1)*(i_points_[j].local_bcoords_A()(2) - i_points_[0].local_bcoords_A()(2)) +
				 i_points_[j].local_bcoords_A()(1)*(i_points_[0].local_bcoords_A()(2) - i_points_[j-1].local_bcoords_A()(2)));
	}
	return fabs(subtotal/2);
};

ostream& operator<<(ostream& os, const IntersectionPolygon& polygon)
{
    for(unsigned int i = 0; i < polygon.i_points_.size(); i++)
        os << polygon.i_points_[i];
    
    return os;
}


} // END NAMESPACE
