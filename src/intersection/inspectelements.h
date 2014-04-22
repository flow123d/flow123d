/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#ifndef INSPECTELEMENTS_H_
#define INSPECTELEMENTS_H_

#include "system.hh"
#include "computeintersection.h"

using namespace computeintersection;
using namespace std;

class InspectElements {

	std::vector<IntersectionLocal *> all_intersections;

	IntersectionLocal temporary_intersection;

	Simplex<1> abscissa;
	Simplex<2> triangle;
	Simplex<3> tetraheadron;

	Mesh *mesh;

	ComputeIntersection<Simplex<1>,Simplex<3> > CI13;
	ComputeIntersection<Simplex<2>,Simplex<3> > CI23;

public:
	InspectElements();
	InspectElements(Mesh *_mesh);
	virtual ~InspectElements();

	void ComputeIntersections23();
	void ComputeIntersections13();

	void UpdateAbscissa(const ElementFullIter &el);
	void UpdateTriangle(const ElementFullIter &el);
	void UpdateTetraheadron(const ElementFullIter &el);
};

#endif /* INSPECTELEMENTS_H_ */
