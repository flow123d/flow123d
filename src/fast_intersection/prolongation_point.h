/*
 * ProlongationPoint.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include <vector>

namespace fast_1_3{

class ProlongationPoint {

	// Index of 1D element
	unsigned int element_1D_idx;

	// Index of 3D element
	unsigned int element_3D_idx;

	// Side index of 3D element
	unsigned int local_side_3D;

	// Local coordinates of side of 3D element
	std::vector<double> coordinates_3D;
	//double alfa;
	//double beta;

	// Local coordinate of 1D element
	double theta;

	// Orientation of 1D element true = right orientation -> theta = theta ---- false = wrong orientation, theta = 1 - theta
	bool orientation;

public:
	ProlongationPoint(unsigned int el1D, unsigned int el3D, unsigned int side3D, std::vector<double> &coords_3D, double t, bool orientation);
	~ProlongationPoint();

	inline unsigned int idx_elm1D(){return element_1D_idx;};
	inline unsigned int idx_elm3D(){return element_3D_idx;};
	inline unsigned int idx_side3D(){return local_side_3D;};
	inline std::vector<double> local_coords_3D(){return coordinates_3D;};
	inline std::vector<double> &local_coords_3D_ref(){return coordinates_3D;};
	inline double local_coords_1D(){return theta;};
	inline double local_coords_1D_if(){
		if(!orientation){ return 1 - theta;}
		return theta;
	}
	inline bool getOrientation(){return orientation;};
};

} // namespace fast_1_3 close
