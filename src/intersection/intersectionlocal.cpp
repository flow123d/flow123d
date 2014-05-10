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


//void IntersectionLocal::add_local_coord(const std::vector<double> &coordin1, const double &coordin2){
	/*for(unsigned int i = 0; i < coordin1.size(); i++){
	std::printf("Souřadnice - 3D: %f \n", coordin1[i]);
	}
	std::printf("Souřadnice - 1D: %f \n", coordin2);
	 */
	//i_points.push_back(new IntersectionPoint(coordin1, coordin2));
//};


} // namespace computeintersection close
