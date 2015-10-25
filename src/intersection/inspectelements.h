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
#include "prolongationpoint.h"
#include "prolongationline.h"
#include "intersectionpolygon.h"



using namespace std;
namespace computeintersection {

/**
* Main class, which takes mesh and you can call method for computing intersections for different dimensions of elements
* It can compute whole polygon area.
* It can create a mesh file with intersections
*/
class InspectElements {

	// For case 2D-3D - list of intersectionpolygon
	std::vector<std::vector<IntersectionPolygon>> intersection_list;
	// For case 1D-3D - list of intersectionline
	std::vector<std::vector<IntersectionLine>> intersection_line_list;
	// Array of flags, which elements are computed
	std::vector<bool> closed_elements;
	std::vector<int> flag_for_3D_elements;


	std::queue<ProlongationLine> prolongation_line_queue_2D;
	std::queue<ProlongationLine> prolongation_line_queue_3D;

	std::queue<ProlongationPoint> prolongation_point_queue;

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

	/**
	 * Every method needs to be implemented for different type of mesh intersection
	 */
	template<unsigned int subdim, unsigned int dim> inline void compute_intersections(){
		cout << "Warning - method compute_intersections "<< subdim<<"D with "<< dim<<"D is not implemented yet" << endl;
	};

	template<unsigned int subdim, unsigned int dim> inline void compute_intersections_init(){
		cout << "Warning - method compute_intersections_init "<< subdim<<"D with "<< dim<<"D is not implemented yet" << endl;
	};


	void computeIntersections2d3d(); // set private
	void computeIntersections2d3dProlongation(const ProlongationLine &pl); // set private
	void computeIntersections2d3dUseProlongationTable(std::vector<unsigned int> &prolongation_table, const ElementFullIter &elm, const ElementFullIter &ele); // set private

	bool intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx); // set private

	void print_mesh_to_file(string name);
	void print_mesh_to_file_1D(string name);

    /** @brief Computes the area of 2d-3d polygonal intersection.
     * @return the polygonal area
     * TODO: probably remove dependence on @class TTriangle in NGH
     */
	double polygonArea();

    /** @brief Computes the length of 1d-3d line intersection.
     * @return the line length
     * TODO: probably remove dependence on @class TAbscissa in NGH
     */
    double line_length();
};

} // END namespace
