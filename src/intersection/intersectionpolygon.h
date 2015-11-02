/*
 * IntersectionPolygon.h
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#ifndef INTERSECTIONPOLYGON_H_
#define INTERSECTIONPOLYGON_H_

#include "intersectionpoint.h"
#include "prolongation.h"
#include "system/system.hh"
#include "mesh/mesh.h"
#include <queue>

using namespace std;
namespace computeintersection{

    /**
     * IntersectionPolygon represents a polygon from computing intersection of a triangle(Simplex<2>) and a tetrahedron(Simplex<3>)
     * Contains array of points (necessary not in order at first)
     * Flag if its "patological" (polygon is a part of side in tetrahedron)
     */
class IntersectionPolygon {

	/// Can be replaced with global epsilon
	static const double epsilon;
	bool is_patological_;

	std::vector<IntersectionPoint<2,3>> i_points;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;

public:
	IntersectionPolygon();
	~IntersectionPolygon();

	IntersectionPolygon(unsigned int elem2D,unsigned int elem3D);

	/**
	 * @param InPoint - other point of polygon inserting into the array
	 */
	void add_ipoint(const IntersectionPoint<2,3> &InPoint);

	inline bool is_patological(){
		return is_patological_;
	};

	inline const IntersectionPoint<2,3> &get_point(const unsigned int index) const
	{
		 return i_points[index];
	}

	inline unsigned int size(){
		return i_points.size();
	}
	inline unsigned int idx_2D(){return element_2D_idx;}
	inline unsigned int idx_3D(){return element_3D_idx;}

	/** method compute local polygon area from barycentric coordinates
	 * @return double computed local area
	 * split the polygon into triangles according to the first point
	 * for every triangle compute area from barycentric coordinates
	 * 				Barycentric coords. 		Local coords.
	 * Point A		[A0,A1,A2]					[A1,A2,0]
	 * Point B	    [B0,B1,B2]					[B1,B2,0]
	 * Point C		[C0,C1,C2]					[C1,C2,0]
	 *
	 *  triangle area: (B-A)x(C-A)/2 => ((B1-A1)(C2-A2) - (B2-A2)(C1-A1))/2
	 *  => (A1(B2-C2) + B1(C2-A2) + C1(A2-B2))/2
	 *
	 *  final polygon area is sum of triangle areas
	 */
	double get_area() const;

	inline void print(){
		for(unsigned int i = 0; i < i_points.size();i++){
			i_points[i].print();
		}
	}

	/**
	 * Tracing generic polygon - if flag "is_patological_" is true - traces polygon by convex hull,
	 * otherwise by optimize method
	 * it fills prolongation table
	 */
	void trace_generic_polygon(std::vector<unsigned int> &prolongation_table);

	void trace_polygon_opt(std::vector<unsigned int> &prolongation_table);

	/// TODO: link for wikipedia
	void trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table);

	double convex_hull_cross(const IntersectionPoint<2,3> &O,
			const IntersectionPoint<2,3> &A,
			const IntersectionPoint<2,3> &B) const;

	/**
	 * @return index of prolongation table for two points, which are on same line of triangle or side of tetrahedron
	 */
	int convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
			const IntersectionPoint<2,3> &B) const;

	/**
	 * @return index of side of tetrahedron - if tetrahedron has polygon on it
	 */
	int side_content_prolongation() const;
};


} // END NAMESPACE
#endif /* INTERSECTIONPOLYGON_H_ */
