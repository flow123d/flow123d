/*
 * inspectelements.h
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 *
 */

#include "system/system.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "fields/field_interpolated_p0.hh"
#include "system/sys_profiler.hh"

#include "fast_intersection/intersection_local.h"
#include "fast_intersection/prolongation_point.h"
#include "bp_fris/intersections.h"
#include <queue>

namespace fast_1_3{

class InspectElements {

public:
	InspectElements();
	/* Constructor takes a pointer mesh and creates vector of bool variables
	 * it calls method calculate_intersections()
	 */
	InspectElements(Mesh* sit);

	/* Main computational method
	 * */
	void calculate_intersections();
	/* Searching method, need to find first Prolongation Point and add it to the queue
	 * */
	bool calculate_prolongation_point(const ElementFullIter &element_1D, const ElementFullIter &element_3D);
	/* Calculating intersections from PPoint
	 * */
	void calculate_from_prolongation_point(ProlongationPoint &point);
	void calculate_intersection_from_1D(unsigned int idx_1D, unsigned int idx_3D, std::vector<double> &interpolated_3D_coords);

	// pomocné výpočetní metody:
	double get_local_coords_1D(SPoint<3> a, SPoint<3> b, SPoint<3> x);
	SPoint<3> local_interpolation(SPoint<3> &point_A,SPoint<3> &point_B,const double &t_A,const double &t_B,const double &t);
	std::vector<double> local_vector_interpolation(const std::vector<double> &point_A,const std::vector<double> &point_B,const double &t_A, const double &t_B, const double &t);
	~InspectElements();

	vector<IntersectionLocal> getIntersections();

	void update_tetrahedron(const ElementFullIter &element_3D);
	void update_abscissa(const ElementFullIter &element_1D, bool orientace);

private:
	// information of all elements if element was inspected
	std::vector<bool> projeti;

	std::vector<IntersectionLocal> all_intersection;

	std::queue<ProlongationPoint> ppoint;

	Mesh* sit;

	// auxiliary data members
	Simplex<3,3> tetrahedron;
	HyperPlane<1,3> abscissa;

};

} // namespace fast_1_3 close
