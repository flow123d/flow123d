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

template<unsigned int N, unsigned int M> class IntersectionPoint;
class IntersectionLine;
class IntersectionPolygon;
    
/**
* Main class, which takes mesh and you can call method for computing intersections for different dimensions of elements
* It can compute whole polygon area.
* It can create a mesh file with intersections
*/
class InspectElements {

    // For case 1D-2D - list of intersection points
    std::vector<std::vector<IntersectionPoint<1,2>>> intersection_point_list;
    
	// For case 2D-3D - list of intersectionpolygon
	std::vector<std::vector<IntersectionPolygon>> intersection_list;
	// For case 1D-3D - list of intersectionline
	std::vector<std::vector<IntersectionLine>> intersection_line_list;
	// Array of flags, which elements are computed
	std::vector<bool> closed_elements;
	std::vector<int> flag_for_3D_elements;

    // Used only for 2d-3d
	std::queue<ProlongationLine> prolongation_line_queue_2D;
	std::queue<ProlongationLine> prolongation_line_queue_3D;

	std::queue<ProlongationPoint> prolongation_point_queue; // Used only for 1d-3d

	int element_2D_index;   // Used only for 2d-3d

	Simplex<1> abscissa;    // Used only for 1d-3d
	Simplex<2> triangle;    // Used only for 2d-3d
	Simplex<3> tetrahedron;

	Mesh *mesh;
	std::vector<BoundingBox> elements_bb;
	BoundingBox mesh_3D_bb;

	void update_abscissa(const ElementFullIter &el);
	void update_triangle(const ElementFullIter &el);
	void update_tetrahedron(const ElementFullIter &el);

    // Used only for 1d-3d
	void prolongate_elements(const IntersectionLine &il, const ElementFullIter &elm, const ElementFullIter &ele);
	void prolongate_1D_element(const ElementFullIter &elm, const ElementFullIter &ele);
	void prolongate(const ProlongationPoint &pp);
    
    // Used only for 2d-3d
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

    const std::vector<IntersectionPoint<1,2>> & list_intersection_points(unsigned int ele_idx);
    const std::vector<IntersectionLine> & list_intersection_lines(unsigned int idx_component_1D);
    
	void print_mesh_to_file(std::string name);
	void print_mesh_to_file_1D(std::string name);

    /** @brief Computes the area of 2d-3d polygonal intersections (sum over all polygons).
     * @return the area of intersection polygon
     */
	double polygon_area();

    /** @brief Computes the length of 1d-3d line intersection (sum over all lines).
     * @return the line length
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


inline const std::vector< IntersectionLine >& InspectElements::list_intersection_lines(unsigned int ele_idx)
{
    ASSERT_LESS(ele_idx, intersection_line_list.size());
    return intersection_line_list[ele_idx];
}

inline const vector< IntersectionPoint< 1, 2 > >& InspectElements::list_intersection_points(unsigned int ele_idx)
{
    ASSERT_LESS(ele_idx, intersection_point_list.size());
    return intersection_point_list[ele_idx];
}


} // END namespace


#endif // INSPECT_ELEMENTS_H_