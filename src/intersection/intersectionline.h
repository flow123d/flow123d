/*
 * IntersectionLine.h
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#ifndef INTERSECTIONLINE_H_
#define INTERSECTIONLINE_H_

#include "intersectionpoint.h"
#include "system/system.hh"
#include "mesh/mesh.h"
#include <queue>

using namespace std;
namespace computeintersection{

class IntersectionLine{

	std::vector<IntersectionPoint<1,3>> i_points;
	unsigned int element_1D_idx;
	unsigned int element_3D_idx;

public:

	inline IntersectionLine(){};
	inline IntersectionLine(unsigned int ele_1D, unsigned int ele_3D):element_1D_idx(ele_1D), element_3D_idx(ele_3D){};
	inline ~IntersectionLine(){};

	inline std::vector<IntersectionPoint<1,3>> &get_points(){
		return i_points;
	}

	inline const std::vector<IntersectionPoint<1,3>> &get_points() const{
		return i_points;
	}

	inline void trace_line(){

		if(i_points.size() > 1){
			unsigned int j = (i_points[0].get_side2() + i_points[0].get_orientation())%2;

			if(j == 0){
				std::vector<IntersectionPoint<1,3>> new_points(2);
				new_points[0] = i_points[1];
				new_points[1] = i_points[0];
				i_points = new_points;
			}
		}

	};

	inline const IntersectionPoint<1,3> &operator[](unsigned int index) const{
		return i_points[index];
	};

	inline const unsigned int size() const{
		return i_points.size();
	};

	inline void add_points(const std::vector<IntersectionPoint<1,3>> &points){
		i_points = points;
	};

	inline unsigned int get_elm1D_idx() const{
		return element_1D_idx;
	};

	inline unsigned int get_elm3D_idx() const{
		return element_3D_idx;
	};

	inline const IntersectionPoint<1,3> &get_point(const unsigned int index) const
	{
		 return i_points[index];
	};

};

} // END NAMESPACE
#endif /* INTERSECTIONLINE_H_ */
