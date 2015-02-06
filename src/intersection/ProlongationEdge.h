/*
 * ProlongationEdge.h
 *
 *  Created on: 17.11.2014
 *      Author: viktor
 */

#include "intersectionpoint.h"
#include "plucker.h"
#include "system/system.hh"

using namespace std;
namespace computeintersection{

class ProlongationEdge {

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;

	bool is_patological;
	std::vector<IntersectionPoint<2,3>> i_points; //vektor ukazatelu na dvojice lokal. souradnic

	// Predani pluckerovych souradnic
	std::vector<Plucker *> p_triangle_coordinates;
	std::vector<Plucker *> p_tetrahedron_coordinates;

	// Predani pluckerovych soucinu

	// Predani objektu computerintersections



public:
	ProlongationEdge();
	virtual ~ProlongationEdge();
};

}
