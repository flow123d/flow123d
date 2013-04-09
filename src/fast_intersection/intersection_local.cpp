/*
 * IntersectionLocal.cpp
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include "intersection_local.h"
#include "prolongation_point.h"

namespace fast_1_3{


IntersectionLocal::IntersectionLocal(unsigned int elem1D,unsigned int elem3D):element_1D_idx(elem1D),element_3D_idx(elem3D){

};

IntersectionLocal::~IntersectionLocal(){};

void IntersectionLocal::add_local_coord(const std::vector<double> &coordin1, const double &coordin2){
	i_points.push_back(new IntersectionPoint(coordin1, coordin2));
};


} // namespace fast_1_3 close
