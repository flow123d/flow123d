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

	/* Possibly replace IntersectionLocal by index of ProlongationLine
	 * in the prolongation_line_queue, which should rather be std::deque
	 * to support both queue and random access operations.
	 */
	std::vector<std::vector<IntersectionLocal>> intersection_list;
	std::vector<std::vector<IntersectionLine>> intersection_line_list;
	std::vector<bool> closed_elements;
	std::vector<int> flag_for_3D_elements;


	std::queue<ProlongationLine> prolongation_line_queue_2D;
	std::queue<ProlongationLine> prolongation_line_queue_3D;

	std::queue<ProlongationPoint> prolongation_point_queue;

	//IntersectionLocal temporary_intersection;

	int element_2D_index;

	Simplex<1> abscissa;
	Simplex<2> triangle;
	Simplex<3> tetrahedron;

	Mesh *mesh;
	std::vector<BoundingBox> elements_bb;
	BoundingBox mesh_3D_bb;

	void update_abscissa(const ElementFullIter &el);
	void update_triangle(const ElementFullIter &el);
	void update_tetrahedron(const ElementFullIter &el);

	void prolongate_elements(const IntersectionLine &il, const ElementFullIter &elm, const ElementFullIter &ele);
	void prolongate_1D_element(const ElementFullIter &elm, const ElementFullIter &ele);
	void prolongate(const ProlongationPoint &pp);

public:

	InspectElements(Mesh *_mesh);
	~InspectElements();

	template<unsigned int subdim, unsigned int dim> inline void compute_intersections(){
		cout << "Warning - method compute_intersections "<< subdim<<"D with "<< dim<<"D is not implemented yet" << endl;
	};

	template<unsigned int subdim, unsigned int dim> inline void compute_intersections_init(){
		cout << "Warning - method compute_intersections_init "<< subdim<<"D with "<< dim<<"D is not implemented yet" << endl;
	};


	void computeIntersections2d3d();
	void computeIntersections2d3dProlongation(const ProlongationLine &pl);
	void computeIntersections2d3dUseProlongationTable(std::vector<unsigned int> &prolongation_table, const ElementFullIter &elm, const ElementFullIter &ele);

	bool intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx);


	void ComputeIntersections23();
	void ComputeIntersections13();

	void print_mesh_to_file(string name);
	void print_mesh_to_file_1D(string name);

	inline double polygonArea(){
		double subtotal = 0.0;

		for(unsigned int i = 0; i < intersection_list.size(); i++){
			for(unsigned int j = 0; j < intersection_list[i].size();j++){
				Element efi = *mesh->element(intersection_list[i][j].idx_2D());
				 TTriangle t2d(efi);
				 double t2dArea = t2d.GetArea();
				 double localArea = intersection_list[i][j].getArea();//il.getArea();
				 subtotal += 2*localArea*t2dArea;
			}
		}
		return subtotal;
	};

};

} // END namespace
