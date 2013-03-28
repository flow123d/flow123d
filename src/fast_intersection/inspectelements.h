/*
 * inspectelements.h
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 *
 */


#include <gtest/gtest.h>
#include "system/system.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "fields/field_interpolated_p0.hh"
#include "system/sys_profiler.hh"

#include "fast_intersection/IntersectionLocal.h"
#include "fast_intersection/ProlongationPoint.h"
#include <queue>


class InspectElements {

public:
	InspectElements();

	InspectElements(Mesh* sit_);

	void calculate_intersections();

	~InspectElements();

	vector<Intersection_Local> getIntersections();

private:
	// information of all elements if element was inspected
	std::vector<bool> projeti;
	// vector of all intersections
	std::vector<Intersection_Local> all_intersection;

	std::queue<ProlongationPoint> ppoint;

	Mesh* sit;


};
