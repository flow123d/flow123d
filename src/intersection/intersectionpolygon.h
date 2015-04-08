/*
 * IntersectionPolygon.h
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#ifndef INTERSECTIONPOLYGON_H_
#define INTERSECTIONPOLYGON_H_

#include "intersectionpoint.h"
#include "prolongationline.h"
#include "system/system.hh"
#include "mesh/mesh.h"
#include <queue>

using namespace std;
namespace computeintersection{

class IntersectionPolygon {

	static const double epsilon;
	bool is_patological_;

	std::vector<IntersectionPoint<2,3>> i_points;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;

public:
	IntersectionPolygon();
	~IntersectionPolygon();

	IntersectionPolygon(unsigned int elem2D,unsigned int elem3D);


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

	/**
	 * method compute local polygon area from barycentric coordinates
	 * @return double computed local area
	 */
	double get_area() const;

	inline void print(){
		for(unsigned int i = 0; i < i_points.size();i++){
			i_points[i].print();
		}
	}

	void trace_generic_polygon(std::vector<unsigned int> &prolongation_table);

	void trace_polygon_opt(std::vector<unsigned int> &prolongation_table);

	void trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table);

	double convex_hull_cross(const IntersectionPoint<2,3> &O,
			const IntersectionPoint<2,3> &A,
			const IntersectionPoint<2,3> &B) const;

	int convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
			const IntersectionPoint<2,3> &B) const;

	int side_content_prolongation() const;
};


} // END NAMESPACE
#endif /* INTERSECTIONPOLYGON_H_ */
