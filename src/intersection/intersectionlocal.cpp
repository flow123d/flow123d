/*
 * IntersectionLocal.cpp
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include "intersectionlocal.h"

namespace computeintersection{

/*
IntersectionLocal::IntersectionLocal(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){
	id = 0;
};
*/
IntersectionLocal::IntersectionLocal(){};

IntersectionLocal::IntersectionLocal(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){

	for(unsigned int i = 0; i < 7;i++){
		for(unsigned int j = 0; j < 3;j++){
			tracing_table(i,j) = -1;
		}
	}
};

IntersectionLocal::~IntersectionLocal(){};

void IntersectionLocal::addIP(IntersectionPoint<2,3> InPoint){
	i_points.push_back(InPoint);
};

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
double IntersectionLocal::getArea(){
	double subtotal = 0.0;
	for(unsigned int j = 2; j < i_points.size();j++){
		//xprintf(Msg, "volani %d %d\n",j, i_points.size());
		subtotal += i_points[0].getLocalCoords1()(1)*(i_points[j-1].getLocalCoords1()(2) - i_points[j].getLocalCoords1()(2)) +
				 i_points[j-1].getLocalCoords1()(1)*(i_points[j].getLocalCoords1()(2) - i_points[0].getLocalCoords1()(2)) +
				 i_points[j].getLocalCoords1()(1)*(i_points[0].getLocalCoords1()(2) - i_points[j-1].getLocalCoords1()(2));
	}
	return subtotal/2;
};


//void IntersectionLocal::add_local_coord(const std::vector<double> &coordin1, const double &coordin2){
	/*for(unsigned int i = 0; i < coordin1.size(); i++){
	std::printf("Souřadnice - 3D: %f \n", coordin1[i]);
	}
	std::printf("Souřadnice - 1D: %f \n", coordin2);
	 */
	//i_points.push_back(new IntersectionPoint(coordin1, coordin2));
//};


} // namespace computeintersection close
