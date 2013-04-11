/*
 * ProlongationPoint.cpp
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include "prolongation_point.h"

namespace fast_1_3{

ProlongationPoint::ProlongationPoint(unsigned int el1D, unsigned int el3D, unsigned int side3D, std::vector<double> &coords_3D,
		double t, bool orientation):element_1D_idx(el1D), element_3D_idx(el3D), local_side_3D(side3D), coordinates_3D(coords_3D), theta(t), orientation(orientation){};







ProlongationPoint::~ProlongationPoint(){};
} // namespace fast_1_3 close
