/*
 * IntersectionLocal.cpp
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include <algorithm>
#include "intersectionlocal.h"

namespace computeintersection{

const double IntersectionLocal::epsilon = 0.000000001;//128*2048*numeric_limits<double>::epsilon();


IntersectionLocal::IntersectionLocal(){};

IntersectionLocal::IntersectionLocal(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){

	is_patological = false;

	/*for(unsigned int i = 0; i < 7;i++){
		for(unsigned int j = 0; j < 4;j++){
			tracing_table(i,j) = -1;
		}
	}*/

	i_points.reserve(4);
};

IntersectionLocal::~IntersectionLocal(){};

void IntersectionLocal::addIP(const IntersectionPoint<2,3> &InPoint){
	if(InPoint.isPatological()){
		is_patological = true;
	}
	i_points.push_back(InPoint);
};

void IntersectionLocal::traceGenericPolygon(std::vector<unsigned int> &prolongation_table){

	if(!is_patological){

		if(element_2D_idx == 5775 && element_3D_idx == 172349){
			//xprintf(Msg, "patrame dal \n");
		}

		trace_polygon_opt(prolongation_table);

		if(element_2D_idx == 5775 && element_3D_idx == 172349){
					//xprintf(Msg, "dopatrano dal \n");
				}
		//xprintf(Msg,"Tracing opt polygon(%d)\n", i_points.size());
		//std::vector<std::pair<unsigned int, unsigned int>> pp;
		//tracePolygonOpt(pp);


	}else{
		if(element_2D_idx == 5775 && element_3D_idx == 172349){
							//xprintf(Msg, "co to kurva \n");
						}
		//xprintf(Msg,"Tracing traceConvexHull polygon(%d)", i_points.size());
		trace_polygon_convex_hull(prolongation_table);
	}
};

void IntersectionLocal::trace_polygon_opt(std::vector<unsigned int> &prolongation_table){

	if(element_2D_idx == 5775 && element_3D_idx == 172349){
		//xprintf(Msg, "trace_polygon_opt %d %d\n", i_points.size(), is_patological);
	}

	if(is_patological || i_points.size() < 2){
		return;
	}

	arma::Mat<int>::fixed<7,2> trace_table;
	for(unsigned int i = 0; i < 7;i++){
		for(unsigned int j = 0; j < 2;j++){
			trace_table(i,j) = -1;
		}
	}

	std::vector<IntersectionPoint<2,3>> new_points;
	new_points.reserve(i_points.size());
	prolongation_table.reserve(i_points.size());
	if(element_2D_idx == 5775 && element_3D_idx == 172349){
			//xprintf(Msg, "trace_polygon_opt2\n");
		}
	/*
	 * Point orientation
	 * S0 + 1 => IN , S0 + 0 => OUT
	 * S1 + 0 => IN
	 * S2 + 1 => IN
	 * S3 + 0 => IN
	 *
	 * H1 has opossite orientation
	 * */

	for(unsigned int i = 0; i < i_points.size();i++){

		if(element_2D_idx == 5775 && element_3D_idx == 172349){
			//i_points[i].print();
				}

			// TYP bodu H-S nebo H-H
			if(i_points[i].getSide1() != -1){
				// TYP bodu H-S
				if(!i_points[i].isVertex()){

					unsigned int j = (i_points[i].getSide2() + i_points[i].getOrientation()+ i_points[i].getSide1())%2;


					if(j == 1){
						trace_table(i_points[i].getSide2(),0) = i_points[i].getSide1()+4;
						trace_table(i_points[i].getSide2(),1) = i;
					}else{
						// bod je vstupní na hraně a vchází do čtyřstěnu
						trace_table(4+i_points[i].getSide1(),0) = i_points[i].getSide2();
						trace_table(4+i_points[i].getSide1(),1) = i;
					}



				}else{
					// TYP bodu H-H
					unsigned int vertex_index = (i_points[i].get_local_coords1()[0] == 1 ? 0 : (i_points[i].get_local_coords1()[1] == 1 ? 1 : 2));
					unsigned int triangle_side_in = RefSimplex<2>::line_sides[vertex_index][1];
					unsigned int triangle_side_out = RefSimplex<2>::line_sides[vertex_index][0];
					if(trace_table(4+triangle_side_in,1) == -1){
						trace_table(4+triangle_side_in,0) = triangle_side_out+4;
						trace_table(4+triangle_side_in,1) = i;
					}
				}
			}else{
				// TYP bodu S-S
				unsigned int tetrahedron_line = i_points[i].getSide2();
				unsigned int tetrahedron_side_in = RefSimplex<3>::line_sides[tetrahedron_line][i_points[i].getOrientation()];
				unsigned int tetrahedron_side_out = RefSimplex<3>::line_sides[tetrahedron_line][1 - i_points[i].getOrientation()];
				trace_table(tetrahedron_side_in,0) = tetrahedron_side_out;
				trace_table(tetrahedron_side_in,1) = i;
			}
	}

	if(element_2D_idx == 5775 && element_3D_idx == 172349){
				//xprintf(Msg, "trace_polygon_opt3\n");
				//trace_table.print();
			}

	//trace_table.print();

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

	prolongation_table.push_back((unsigned int)trace_table(first_row_index,0));

	while(first_row_index != next_row){
		new_points.push_back(i_points[(unsigned int)trace_table(next_row,1)]);
		//prolongation_table.push_back((unsigned int)next_row);
		prolongation_table.push_back((unsigned int)trace_table(next_row,0));
		next_row = trace_table(next_row,0);
	}
	//prolongation_table.push_back((unsigned int)next_row);

	i_points = new_points;

};

void IntersectionLocal::trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table){


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
	if((int)i_points.size() != n){
		xprintf(Msg, "Velikost před a po (%d, %d)\n", n, i_points.size());

	}

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

double IntersectionLocal::convex_hull_cross(const IntersectionPoint<2,3> &O,
		const IntersectionPoint<2,3> &A,
		const IntersectionPoint<2,3> &B) const{
	return ((A.get_local_coords1()[1]-O.get_local_coords1()[1])*(B.get_local_coords1()[2]-O.get_local_coords1()[2])
			-(A.get_local_coords1()[2]-O.get_local_coords1()[2])*(B.get_local_coords1()[1]-O.get_local_coords1()[1]));
}

int IntersectionLocal::convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
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
	cout << A.get_local_coords1()[0] << ","
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
	xprintf(Msg, "Mozny problem\n");
	return -1;
}

int IntersectionLocal::side_content_prolongation() const{

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

/* split the polygon into triangles according to the first point
 * for every triangle compute area from barycentric coordinates
 * 				Barycentric coords. 		Local coords.
 * Point A		[A0,A1,A2]					[A1,A2,0]
 * Point B	    [B0,B1,B2]					[B1,B2,0]
 * Point C		[C0,C1,C2]					[C1,C2,0]
 *
 *  triangle area: (B-A)x(C-A)/2 => ((B1-A1)(C2-A2) - (B2-A2)(C1-A1))/2
 *  => (A1(B2-C2) + B1(C2-A2) + C1(A2-B2))/2
 *
 *  final polygon area is sum of triangle areas
 */
double IntersectionLocal::getArea() const{
	double subtotal = 0.0;
	for(unsigned int j = 2; j < i_points.size();j++){
		//xprintf(Msg, "volani %d %d\n",j, i_points.size());
		subtotal += fabs(i_points[0].get_local_coords1()(1)*(i_points[j-1].get_local_coords1()(2) - i_points[j].get_local_coords1()(2)) +
				 i_points[j-1].get_local_coords1()(1)*(i_points[j].get_local_coords1()(2) - i_points[0].get_local_coords1()(2)) +
				 i_points[j].get_local_coords1()(1)*(i_points[0].get_local_coords1()(2) - i_points[j-1].get_local_coords1()(2)));
	}
	return fabs(subtotal/2);
};


} // namespace computeintersection close
