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

	std::map<unsigned int, std::vector<IntersectionLocal>> intersection_list;
	std::map<unsigned int, bool> closed_elements;

	std::queue<ProlongationLine> prolongation_line_queue;

	IntersectionLocal temporary_intersection;

	Simplex<1> abscissa;
	Simplex<2> triangle;
	Simplex<3> tetrahedron;

	Mesh *mesh;

	ComputeIntersection<Simplex<1>,Simplex<3> > CI13;
	ComputeIntersection<Simplex<2>,Simplex<3> > CI23;

public:
	InspectElements();
	InspectElements(Simplex<2> sim2, Simplex<3> sim3);
	InspectElements(Mesh *_mesh);
	~InspectElements();

	// dopnit implementaci do .cpp
	void computeIntersections2d3d();
	void computeIntersections2d3dInit();
	void computeIntersections2d3dProlongation();


	void ComputeIntersections23();
	void ComputeIntersections13();

	void UpdateAbscissa(const ElementFullIter &el);
	void UpdateTriangle(const ElementFullIter &el);
	void UpdateTetrahedron(const ElementFullIter &el);

	void print(unsigned int vyber);
	inline double polygonArea(){
		double subtotal = 0.0;
		for(unsigned int i = 0; i < all_intersections.size();i++){
			 Element efi = *mesh->element(all_intersections[i].idx_2D());
			 TTriangle t2d(efi);
			 double t2dArea = t2d.GetArea();
			 double localArea = all_intersections[i].getArea();//il.getArea();
			 subtotal += 2*localArea*t2dArea;
			 xprintf(Msg,"Subtotal: %f\n",2*localArea*t2dArea);
		}
		return subtotal;
	}

};

} // END namespace
