/*
 * inspectelements.h
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */
#ifndef INSPECT_ELEMENTS_H_
#define INSPECT_ELEMENTS_H_

#include "system/system.hh"

#include "mesh/bounding_box.hh"
#include "mesh/mesh_types.hh"

#include "simplex.h"

#include <queue>

class Mesh; // forward declare

namespace computeintersection {

// forward declare
struct ProlongationLine;
struct ProlongationPoint;
class IntersectionLine;
class IntersectionPolygon;
    
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
    
    void computeIntersections2d3d();
    void computeIntersections2d3dProlongation(const ProlongationLine &pl);
    void computeIntersections2d3dUseProlongationTable(std::vector<unsigned int> &prolongation_table, 
                                                      const ElementFullIter &elm, 
                                                      const ElementFullIter &ele);

    bool intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx);

public:

	InspectElements(Mesh *_mesh);
	~InspectElements();

	/**
	 * Every method needs to be implemented for different type of mesh intersection
	 */
	template<unsigned int subdim, unsigned int dim> 
	void compute_intersections();

	template<unsigned int subdim, unsigned int dim>
	void compute_intersections_init();

	void print_mesh_to_file(std::string name);
	void print_mesh_to_file_1D(std::string name);

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





template<unsigned int subdim, unsigned int dim> 
void InspectElements::compute_intersections()
{
        ASSERT(0, "Method not implemented for given dim and subdim.");
}

template<unsigned int subdim, unsigned int dim> 
void InspectElements::compute_intersections_init()
{
        ASSERT(0, "Method not implemented for given dim and subdim.");
}

} // END namespace


#endif // INSPECT_ELEMENTS_H_