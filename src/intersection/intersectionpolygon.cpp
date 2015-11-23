/*
 * IntersectionPolygon.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "intersectionpolygon.h"
#include "mesh/ref_element.hh"
#include <algorithm>

namespace computeintersection{

const double IntersectionPolygon::epsilon = 0.000000001;

IntersectionPolygon::IntersectionPolygon() {};

IntersectionPolygon::~IntersectionPolygon() {
	// TODO Auto-generated destructor stub
}

IntersectionPolygon::IntersectionPolygon(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){

	is_patological_ = false;
	i_points.reserve(4);
};

void IntersectionPolygon::add_ipoint(const IntersectionPoint<2,3> &InPoint){
	if(InPoint.is_patological()){
		is_patological_ = true;
	}
	i_points.push_back(InPoint);
};

void IntersectionPolygon::trace_generic_polygon(std::vector<unsigned int> &prolongation_table){

	if(!is_patological_){
		trace_polygon_opt(prolongation_table);
	}else{
		trace_polygon_convex_hull(prolongation_table);
	}
};

void IntersectionPolygon::trace_polygon_opt(std::vector<unsigned int> &prolongation_table){

    ASSERT(!is_patological_, "Cannot call this polygonal tracing in pathologic case.");
    
    // avoid tracing (none is needed) if the intersection is just single point
	if(i_points.size() < 2) return;

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
	new_points.reserve(i_points.size());
	prolongation_table.reserve(i_points.size());

//     for(unsigned int i = 0; i < i_points.size();i++)
//         i_points[i].print();

    // go through all intersection points (vertices of polygon)
	for(unsigned int i = 0; i < i_points.size();i++){
            IntersectionPoint<2,3> ip = i_points[i];
			if(ip.get_side1() != -1 ){ 
                // if the edge of triangle is set, it must be
                // intersection S-E or E-S or E-E
                
				if(!ip.is_vertex()){   // IP is not a vertex of triangle -> cannot be E-E
                    DBGMSG("Intersection type E%d - S%d, ori %d.\n",ip.get_side1(),ip.get_side2(),ip.get_orientation());
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
					unsigned int j = (RefElement<3>::normal_orientation(ip.get_side2()) + 
                                      RefElement<2>::normal_orientation(ip.get_side1()) +
                                      ip.get_orientation()
                                     ) % 2;

					if(j == 1){
						// bod je vstupní na stěně a pokračuje na hranu trojúhelníku
                        unsigned int row = ip.get_side2();
						trace_table(row,0) = ip.get_side1()+4;
						trace_table(row,1) = i;
                        DBGMSG("S-E: row: %d, to edge: %d, ip: %d \n",row, ip.get_side1()+4, i);
					}else{
						// bod je vstupní na hraně a vchází do čtyřstěnu
                        unsigned int row = 4+ip.get_side1();
						trace_table(row,0) = ip.get_side2();
						trace_table(row,1) = i;
                        DBGMSG("E-S: row: %d, to side: %d, ip: %d \n",row, ip.get_side2(), i);
					}


				}else{  // type E-E
                    DBGMSG("Intersection type E-E.\n");
                    // IP is the vertex of triangle,
                    // the directions of triangle edges determines the directions of polygon edges
                    //TODO we must provide a function in reference element: node coordinates -> node index
                    //TODO even better: we should save vertex_index when computed
// 					unsigned int vertex_index = (ip.get_local_coords1()[0] == 1 ? 0 : (ip.get_local_coords1()[1] == 1 ? 1 : 2));
                    // using ref element, but we still want to get rid of this...
                    unsigned int vertex_index;
                    RefElement<2>::vertex_index(ip.get_local_coords1(),vertex_index);
                    DBGMSG("E-E: vertex_index: %d \n",vertex_index);
                    // <2>::line_sides[vertex_index][0 = line index IN, or 1 = line index OUT]
					unsigned int triangle_side_in = RefElement<2>::line_sides[vertex_index][0];
					unsigned int triangle_side_out = RefElement<2>::line_sides[vertex_index][1];
                    unsigned int row = 4+triangle_side_in;
                    
                    if(trace_table(row,1) == -1)    // this ignores duplicit IP
                    {
                        trace_table(row,0) = triangle_side_out+4;
                        trace_table(row,1) = i;
                        DBGMSG("E-E: row: %d, to edge: %d, ip: %d \n",row, triangle_side_out+4, i);
                    }
				}
			}else{  // type S-S, IP is on the edge between two sides
                DBGMSG("Intersection type S-S, ori %d.\n", ip.get_orientation());
				/** here the side2 contains number of edge of the tetrahedron where the IP lies
                 * sides: let edge be oriented up and let us see the tetrahedron from outside ->
                 *  then the right side [0] is in and left side [1] is out
                 *  - this is determined by IP orientation
                 */
                unsigned int tetrahedron_line = ip.get_side2();
                unsigned int tetrahedron_side_in = RefElement<3>::line_sides[tetrahedron_line][1-ip.get_orientation()];
                unsigned int tetrahedron_side_out = RefElement<3>::line_sides[tetrahedron_line][ip.get_orientation()];
                
				trace_table(tetrahedron_side_in,0) = tetrahedron_side_out;
				trace_table(tetrahedron_side_in,1) = i;
                DBGMSG("S-S: row: %d, to side: %d, ip: %d \n",tetrahedron_side_in, tetrahedron_side_out, i);
			}
	}
    trace_table.print();
	// Procházení trasovací tabulky a přeuspořádání bodů do nového pole
	int first_row_index = -1;
	int next_row = -1;
	for(unsigned int i = 0; i < 7;i++){
		if(first_row_index != -1){
			if(trace_table(i,0) == first_row_index){
				new_points.push_back(i_points[(unsigned int)trace_table(i,1)]);
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
        ASSERT_LESS(i_ip_orig, i_points.size());
		new_points.push_back(i_points[i_ip_orig]);
		prolongation_table.push_back((unsigned int)trace_table(next_row,0));
		next_row = trace_table(next_row,0);
	}

	i_points = new_points;

};

void IntersectionPolygon::trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table){

	if(i_points.size() <= 1){
			return;
	}

	std::sort(i_points.begin(),i_points.end());

	// Odstranění duplicit
	// Odstranit nejen stejne, ale o epsilon podobne
	for(unsigned int i = 0; i < i_points.size()-1; i++){
		if((fabs(i_points[i].get_local_coords1()[0] - i_points[i+1].get_local_coords1()[0]) < epsilon) &&
		   (fabs(i_points[i].get_local_coords1()[1] - i_points[i+1].get_local_coords1()[1]) < epsilon)){
			i_points.erase(i_points.begin()+i);
		}
	}

	if(i_points.size() > 1 && fabs(i_points[0].get_local_coords1()[0] - i_points[i_points.size()-1].get_local_coords1()[0]) < epsilon &&
			fabs(i_points[0].get_local_coords1()[1] - i_points[i_points.size()-1].get_local_coords1()[1]) < epsilon){
		i_points.erase(i_points.end());
	}


	int n = i_points.size(), k = 0;
	std::vector<IntersectionPoint<2,3>> H(2*n);

	for(int i = 0; i < n; ++i){
		while(k >= 2 && convex_hull_cross(H[k-2], H[k-1], i_points[i]) <= 0) k--;
		H[k++] = i_points[i];
	}

	for(int i = n-2, t = k+1;i>=0;i--){
		while(k >= t && (convex_hull_cross(H[k-2], H[k-1], i_points[i])) <= 0) k--;
		H[k++] = i_points[i];
	}

	H.resize(k-1);
	i_points = H;


	// Filling prolongation table
	if(i_points.size() > 1){

		// Prolongation - finding same zero bary coord/coords except
		// side with content intersect
		unsigned int size = i_points.size() == 2 ? 1 : i_points.size();
		int forbidden_side = side_content_prolongation();

		for(unsigned int i = 0; i < size;i++){
			unsigned int ii = (i+1)%i_points.size();

			for(unsigned int j = 0; j < 3;j++){
				if(fabs((double)i_points[i].get_local_coords1()[j]) < epsilon && fabs((double)i_points[ii].get_local_coords1()[j]) < epsilon){
					prolongation_table.push_back(6-j);
				}
			}
			for(unsigned int k = 0; k < 4;k++){
				if((int)k != forbidden_side && fabs((double)i_points[i].get_local_coords2()[k]) < epsilon && fabs((double)i_points[ii].get_local_coords2()[k]) < epsilon){
					prolongation_table.push_back(3-k);
				}
			}
		}
	}

};

double IntersectionPolygon::convex_hull_cross(const IntersectionPoint<2,3> &O,
		const IntersectionPoint<2,3> &A,
		const IntersectionPoint<2,3> &B) const{
	return ((A.get_local_coords1()[1]-O.get_local_coords1()[1])*(B.get_local_coords1()[2]-O.get_local_coords1()[2])
			-(A.get_local_coords1()[2]-O.get_local_coords1()[2])*(B.get_local_coords1()[1]-O.get_local_coords1()[1]));
}

int IntersectionPolygon::convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
    		const IntersectionPoint<2,3> &B) const{

	// Finding same zero 2D local coord

	if(fabs((double)A.get_local_coords1()[0]) < epsilon && fabs((double)B.get_local_coords1()[0]) < epsilon){
		return 6;
	}else if(fabs((double)A.get_local_coords1()[1]) < epsilon && fabs((double)B.get_local_coords1()[1]) < epsilon){
		return 5;
	}else if(fabs((double)A.get_local_coords1()[2]) < epsilon && fabs((double)B.get_local_coords1()[2]) < epsilon){
		return 4;
	}else if(fabs((double)A.get_local_coords2()[0]) < epsilon && fabs((double)B.get_local_coords2()[0]) < epsilon){
		return 3;
	}else if(fabs((double)A.get_local_coords2()[1]) < epsilon && fabs((double)B.get_local_coords2()[1]) < epsilon){
		return 2;
	}else if(fabs((double)A.get_local_coords2()[2]) < epsilon && fabs((double)B.get_local_coords2()[2]) < epsilon){
		return 1;
	}else if(fabs((double)A.get_local_coords2()[3]) < epsilon && fabs((double)B.get_local_coords2()[3]) < epsilon){
		return 0;
	}
	// Todo: Nahradit Assertem
	/*cout << A.get_local_coords1()[0] << ","
			<< A.get_local_coords1()[1] << ","
			<< A.get_local_coords1()[2] << ","
			<< A.get_local_coords2()[0] << ","
			<< A.get_local_coords2()[1] << ","
			<< A.get_local_coords2()[2] << ","
			<< A.get_local_coords2()[3] << "," << endl;


	cout << B.get_local_coords1()[0] << ","
			<< B.get_local_coords1()[1] << ","
			<< B.get_local_coords1()[2] << ","
			<< B.get_local_coords2()[0] << ","
			<< B.get_local_coords2()[1] << ","
			<< B.get_local_coords2()[2] << ","
			<< B.get_local_coords2()[3] << ",";
	xprintf(Msg, "Mozny problem\n");*/
	return -1;
}

int IntersectionPolygon::side_content_prolongation() const{

	vector<unsigned int> side_field(4,0);

	for(unsigned int i = 0; i < i_points.size();i++){

		for(unsigned int j = 0; j < 4;j++){
			if(fabs((double)i_points[i].get_local_coords2()[j]) < epsilon){
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
	for(unsigned int j = 2; j < i_points.size();j++){
		//xprintf(Msg, "volani %d %d\n",j, i_points.size());
		subtotal += fabs(i_points[0].get_local_coords1()(1)*(i_points[j-1].get_local_coords1()(2) - i_points[j].get_local_coords1()(2)) +
				 i_points[j-1].get_local_coords1()(1)*(i_points[j].get_local_coords1()(2) - i_points[0].get_local_coords1()(2)) +
				 i_points[j].get_local_coords1()(1)*(i_points[0].get_local_coords1()(2) - i_points[j-1].get_local_coords1()(2)));
	}
	return fabs(subtotal/2);
};

ostream& operator<<(ostream& os, const IntersectionPolygon& polygon)
{
    for(unsigned int i = 0; i < polygon.i_points.size(); i++)
        os << polygon.i_points[i];
    
    return os;
}


} // END NAMESPACE
