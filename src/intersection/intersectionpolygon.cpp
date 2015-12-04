/*
 * IntersectionPolygon.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "intersectionpolygon.h"
#include "intersectionpoint.h"
#include "mesh/ref_element.hh"
#include <algorithm>

namespace computeintersection{

const double IntersectionPolygon::epsilon = 0.000000001;

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
    
    // avoid tracing (none is needed) if the intersection is just single point
	if(i_points_.size() < 2) return;

	/**
	 * Vytvoříme tabulku 7x2, první 4 řádky pro stěny čtyřstěnu, další 3 pro hrany trojúhelníku
	 * 1. sloupeček označuje index dalšího řádku na který se pokračuje
	 * 2. sloupeček obsahuje index bodu v původním poli
	 */
	arma::Mat<int>::fixed<7,2> trace_table;
	for(unsigned int i = 0; i < 7;i++){
		for(unsigned int j = 0; j < 2;j++){
			trace_table(i,j) = -1;
		}
	}

	std::vector<IntersectionPoint<2,3>> new_points;
	new_points.reserve(i_points_.size());
	prolongation_table.reserve(i_points_.size());

//     for(unsigned int i = 0; i < i_points_.size();i++)
//         i_points_[i].print();

    // go through all intersection points (vertices of polygon)
	for(unsigned int i = 0; i < i_points_.size();i++){
        IntersectionPoint<2,3> ip = i_points_[i];
        
        switch(ip.dim_A())
        {
            case 0: // if the ip is at the vertex of triangle -> intersection E-E
            {
                DBGMSG("Intersection type E-E.\n");
                // IP is the vertex of triangle,
                // the directions of triangle edges determines the directions of polygon edges
                unsigned int vertex_index = ip.side_idx1();
                DBGMSG("E-E: vertex_index: %d %d\n", vertex_index);
                // <2>::lines_nodes[vertex_index][0 = line index IN, or 1 = line index OUT]
                unsigned int triangle_side_in = RefElement<2>::interact<1,0>(vertex_index)[0];
                unsigned int triangle_side_out = RefElement<2>::interact<1,0>(vertex_index)[1];
                unsigned int row = 4+triangle_side_in;
                
                if(trace_table(row,1) == -1)    // this ignores duplicit IP
                {
                    trace_table(row,0) = triangle_side_out+4;
                    trace_table(row,1) = i;
                    DBGMSG("E-E: row: %d, to edge: %d, ip: %d \n",row, triangle_side_out+4, i);
                }
            } break;
            case 1: // if the ip is on the side of triangle -> intersection S-E or E-S
            {
                DBGMSG("Intersection type E%d - S%d, ori %d.\n",ip.side_idx1(),ip.side_idx2(),ip.orientation());
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
                unsigned int j = (RefElement<3>::normal_orientation(ip.side_idx2()) + 
                                  RefElement<2>::normal_orientation(ip.side_idx1()) +
                                  ip.orientation()
                                 ) % 2;

                if(j == 1){
                    // bod je vstupní na stěně a pokračuje na hranu trojúhelníku
                    unsigned int row = ip.side_idx2();
                    trace_table(row,0) = ip.side_idx1()+4;
                    trace_table(row,1) = i;
                    DBGMSG("S-E: row: %d, to edge: %d, ip: %d \n",row, ip.side_idx1()+4, i);
                }else{
                    // bod je vstupní na hraně a vchází do čtyřstěnu
                    unsigned int row = 4+ip.side_idx1();
                    trace_table(row,0) = ip.side_idx2();
                    trace_table(row,1) = i;
                    DBGMSG("E-S: row: %d, to side: %d, ip: %d \n",row, ip.side_idx2(), i);
                }
            } break;
            case 2: // type S-S, IP is on the edge between two sides
            {
                DBGMSG("Intersection type S-S, ori %d.\n", ip.orientation());
                /** here the side2 contains number of edge of the tetrahedron where the IP lies
                 * sides: let edge be oriented up and let us see the tetrahedron from outside ->
                 *  then the right side [0] is in and left side [1] is out
                 *  - this is determined by IP orientation
                 */
                unsigned int tetrahedron_line = ip.side_idx2();
                unsigned int tetrahedron_side_in = RefElement<3>::interact<2,1>(tetrahedron_line)[1-ip.orientation()];
                unsigned int tetrahedron_side_out = RefElement<3>::interact<2,1>(tetrahedron_line)[ip.orientation()];
                
                trace_table(tetrahedron_side_in,0) = tetrahedron_side_out;
                trace_table(tetrahedron_side_in,1) = i;
                DBGMSG("S-S: row: %d, to side: %d, ip: %d \n",tetrahedron_side_in, tetrahedron_side_out, i);
            } break;
        }
        
    }
            
            
//             if(ip.dim_A() < 2) {
//                 // if the ip is on the side or at the vertex of triangle -> intersection S-E or E-S or E-E
//                 
// 			//if(ip.side_idx1() != unset_loc_idx ){ 
//                 // if the edge of triangle is set, it must be
//                 // intersection S-E or E-S or E-E
//                 
// 				if(!ip.is_vertex()){   // IP is not a vertex of triangle -> cannot be E-E
//                     DBGMSG("Intersection type E%d - S%d, ori %d.\n",ip.side_idx1(),ip.side_idx2(),ip.orientation());
//                     /* directions of polygon edges depends on three aspects:
//                      * - normal orientation of the tetrahedron side (OUT...0, IN...1)
//                      * - orientation of the edge of the triangle (according to direction of other edges of triangle or opposite)
//                      * - orientation of the intersection point ... Plucker product (P>0...1, P<0...0)
//                      * 
//                      * we can obtain this XOR table: XOR(XOR(S,P),E)
//                      *  S   P   E   RES
//                      *  0   0   0   0
//                      *  0   0   1   1
//                      *  0   1   0   1
//                      *  0   1   1   0
//                      *  1   0   0   1
//                      *  1   0   1   0
//                      *  1   1   0   0
//                      *  1   1   1   1
//                      * 
//                      * -> this can be expressed as sum(S,P,E)%2
//                      * 
//                      * if RES == 1 -> intersection direction is S-E (tetrahedron Side -> triangle Edge)
//                      * if RES == 0 -> intersection direction is E-S
//                      */
// 					unsigned int j = (RefElement<3>::normal_orientation(ip.side_idx2()) + 
//                                       RefElement<2>::normal_orientation(ip.side_idx1()) +
//                                       ip.orientation()
//                                      ) % 2;
// 
// 					if(j == 1){
// 						// bod je vstupní na stěně a pokračuje na hranu trojúhelníku
//                         unsigned int row = ip.side_idx2();
// 						trace_table(row,0) = ip.side_idx1()+4;
// 						trace_table(row,1) = i;
//                         DBGMSG("S-E: row: %d, to edge: %d, ip: %d \n",row, ip.side_idx1()+4, i);
// 					}else{
// 						// bod je vstupní na hraně a vchází do čtyřstěnu
//                         unsigned int row = 4+ip.side_idx1();
// 						trace_table(row,0) = ip.side_idx2();
// 						trace_table(row,1) = i;
//                         DBGMSG("E-S: row: %d, to side: %d, ip: %d \n",row, ip.side_idx2(), i);
// 					}
// 
// 
// 				}else{  // type E-E
//                     DBGMSG("Intersection type E-E.\n");
//                     // IP is the vertex of triangle,
//                     // the directions of triangle edges determines the directions of polygon edges
//                     unsigned int vertex_index = ip.side_idx1();
//                     DBGMSG("E-E: vertex_index: %d %d\n", vertex_index);
//                     // <2>::lines_nodes[vertex_index][0 = line index IN, or 1 = line index OUT]
// 					unsigned int triangle_side_in = RefElement<2>::interact<1,0>(vertex_index)[0];
// 					unsigned int triangle_side_out = RefElement<2>::interact<1,0>(vertex_index)[1];
//                     unsigned int row = 4+triangle_side_in;
//                     
//                     if(trace_table(row,1) == -1)    // this ignores duplicit IP
//                     {
//                         trace_table(row,0) = triangle_side_out+4;
//                         trace_table(row,1) = i;
//                         DBGMSG("E-E: row: %d, to edge: %d, ip: %d \n",row, triangle_side_out+4, i);
//                     }
// 				}
// 			}else{  // type S-S, IP is on the edge between two sides
//                 DBGMSG("Intersection type S-S, ori %d.\n", ip.orientation());
// 				/** here the side2 contains number of edge of the tetrahedron where the IP lies
//                  * sides: let edge be oriented up and let us see the tetrahedron from outside ->
//                  *  then the right side [0] is in and left side [1] is out
//                  *  - this is determined by IP orientation
//                  */
//                 unsigned int tetrahedron_line = ip.side_idx2();
//                 unsigned int tetrahedron_side_in = RefElement<3>::interact<2,1>(tetrahedron_line)[1-ip.orientation()];//RefElement<3>::line_sides[tetrahedron_line][1-ip.orientation()];
//                 unsigned int tetrahedron_side_out = RefElement<3>::interact<2,1>(tetrahedron_line)[ip.orientation()];
//                 
// 				trace_table(tetrahedron_side_in,0) = tetrahedron_side_out;
// 				trace_table(tetrahedron_side_in,1) = i;
//                 DBGMSG("S-S: row: %d, to side: %d, ip: %d \n",tetrahedron_side_in, tetrahedron_side_out, i);
// 			}
// 	}
    trace_table.print();
	// Procházení trasovací tabulky a přeuspořádání bodů do nového pole
	int first_row_index = -1;
	int next_row = -1;
	for(unsigned int i = 0; i < 7;i++){
		if(first_row_index != -1){
			if(trace_table(i,0) == first_row_index){
				new_points.push_back(i_points_[(unsigned int)trace_table(i,1)]);
				first_row_index = i;
				next_row = trace_table(i,0);
				break;
			}
		}else if(trace_table(i,0) != -1){
			first_row_index = i;
			next_row = trace_table(i,0);
		}
	}

	// Naplnění prodlužovací tabulky
	// Prodlužovací tabulka obsahuje indexy stěn a (indexy+4) hran
	prolongation_table.push_back((unsigned int)trace_table(first_row_index,0));

	while(first_row_index != next_row){
        DBGMSG("next_row = %d\n",next_row);
        unsigned int i_ip_orig = (unsigned int)trace_table(next_row,1);
        ASSERT_LESS(i_ip_orig, i_points_.size());
		new_points.push_back(i_points_[i_ip_orig]);
		prolongation_table.push_back((unsigned int)trace_table(next_row,0));
		next_row = trace_table(next_row,0);
	}

	i_points_ = new_points;

};

void IntersectionPolygon::trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table){

	if(i_points_.size() <= 1){
			return;
	}

	std::sort(i_points_.begin(),i_points_.end());

	// Odstranění duplicit
	// Odstranit nejen stejne, ale o epsilon podobne
	for(unsigned int i = 0; i < i_points_.size()-1; i++){
		if((fabs(i_points_[i].local_coords1()[0] - i_points_[i+1].local_coords1()[0]) < epsilon) &&
		   (fabs(i_points_[i].local_coords1()[1] - i_points_[i+1].local_coords1()[1]) < epsilon)){
			i_points_.erase(i_points_.begin()+i);
		}
	}

	if(i_points_.size() > 1 && fabs(i_points_[0].local_coords1()[0] - i_points_[i_points_.size()-1].local_coords1()[0]) < epsilon &&
			fabs(i_points_[0].local_coords1()[1] - i_points_[i_points_.size()-1].local_coords1()[1]) < epsilon){
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
				if(fabs((double)i_points_[i].local_coords1()[j]) < epsilon && fabs((double)i_points_[ii].local_coords1()[j]) < epsilon){
					prolongation_table.push_back(6-j);
				}
			}
			for(unsigned int k = 0; k < 4;k++){
				if((int)k != forbidden_side && fabs((double)i_points_[i].local_coords2()[k]) < epsilon && fabs((double)i_points_[ii].local_coords2()[k]) < epsilon){
					prolongation_table.push_back(3-k);
				}
			}
		}
	}

};

double IntersectionPolygon::convex_hull_cross(const IntersectionPoint<2,3> &O,
		const IntersectionPoint<2,3> &A,
		const IntersectionPoint<2,3> &B) const{
	return ((A.local_coords1()[1]-O.local_coords1()[1])*(B.local_coords1()[2]-O.local_coords1()[2])
			-(A.local_coords1()[2]-O.local_coords1()[2])*(B.local_coords1()[1]-O.local_coords1()[1]));
}

// int IntersectionPolygon::convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
//     		const IntersectionPoint<2,3> &B) const{
// 
// 	// Finding same zero 2D local coord
// 
// 	if(fabs((double)A.local_coords1()[0]) < epsilon && fabs((double)B.local_coords1()[0]) < epsilon){
// 		return 6;
// 	}else if(fabs((double)A.local_coords1()[1]) < epsilon && fabs((double)B.local_coords1()[1]) < epsilon){
// 		return 5;
// 	}else if(fabs((double)A.local_coords1()[2]) < epsilon && fabs((double)B.local_coords1()[2]) < epsilon){
// 		return 4;
// 	}else if(fabs((double)A.local_coords2()[0]) < epsilon && fabs((double)B.local_coords2()[0]) < epsilon){
// 		return 3;
// 	}else if(fabs((double)A.local_coords2()[1]) < epsilon && fabs((double)B.local_coords2()[1]) < epsilon){
// 		return 2;
// 	}else if(fabs((double)A.local_coords2()[2]) < epsilon && fabs((double)B.local_coords2()[2]) < epsilon){
// 		return 1;
// 	}else if(fabs((double)A.local_coords2()[3]) < epsilon && fabs((double)B.local_coords2()[3]) < epsilon){
// 		return 0;
// 	}
// 	// Todo: Nahradit Assertem
// 	/*cout << A.local_coords1()[0] << ","
// 			<< A.local_coords1()[1] << ","
// 			<< A.local_coords1()[2] << ","
// 			<< A.local_coords2()[0] << ","
// 			<< A.local_coords2()[1] << ","
// 			<< A.local_coords2()[2] << ","
// 			<< A.local_coords2()[3] << "," << endl;
// 
// 
// 	cout << B.local_coords1()[0] << ","
// 			<< B.local_coords1()[1] << ","
// 			<< B.local_coords1()[2] << ","
// 			<< B.local_coords2()[0] << ","
// 			<< B.local_coords2()[1] << ","
// 			<< B.local_coords2()[2] << ","
// 			<< B.local_coords2()[3] << ",";
// 	xprintf(Msg, "Mozny problem\n");*/
// 	return -1;
// }

int IntersectionPolygon::side_content_prolongation() const{

	vector<unsigned int> side_field(4,0);

	for(unsigned int i = 0; i < i_points_.size();i++){

		for(unsigned int j = 0; j < 4;j++){
			if(fabs((double)i_points_[i].local_coords2()[j]) < epsilon){
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
		subtotal += fabs(i_points_[0].local_coords1()(1)*(i_points_[j-1].local_coords1()(2) - i_points_[j].local_coords1()(2)) +
				 i_points_[j-1].local_coords1()(1)*(i_points_[j].local_coords1()(2) - i_points_[0].local_coords1()(2)) +
				 i_points_[j].local_coords1()(1)*(i_points_[0].local_coords1()(2) - i_points_[j-1].local_coords1()(2)));
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
