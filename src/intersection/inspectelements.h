/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "mesh/mesh.h"
#include "mesh/bih_tree.hh"

#include "mesh/bih_tree.hh"
#include "fields/field_interpolated_p0.hh"

#include "computeintersection.h"

using namespace std;
namespace computeintersection {


class InspectElements {

	std::vector<IntersectionLocal> all_intersections;

	IntersectionLocal temporary_intersection;

	Simplex<1> abscissa;
	Simplex<2> triangle;
	Simplex<3> tetrahedron;

	Mesh *mesh;

	ComputeIntersection<Simplex<1>,Simplex<3> > CI13;
	ComputeIntersection<Simplex<2>,Simplex<3> > CI23;

public:
	InspectElements();
	InspectElements(Mesh *_mesh);
	~InspectElements();

	void ComputeIntersections23();
	void ComputeIntersections13();

	void UpdateAbscissa(const ElementFullIter &el);
	void UpdateTriangle(const ElementFullIter &el);
	void UpdateTetrahedron(const ElementFullIter &el);

	void print(unsigned int vyber);
	inline void polygonArea(){
		for(unsigned int i = 0; i < all_intersections.size();i++){
			IntersectionLocal il = all_intersections[i];

		}
	}

};

} // END namespace
