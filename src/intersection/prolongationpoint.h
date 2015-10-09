#ifndef PROLONGATIONPOINT_H_
#define PROLONGATIONPOINT_H_

using namespace std;
namespace computeintersection{

/**
 * Simple class defines indices of elements for later processing of computing intersection 1D-3D
 * TODO:
 * - replace with struct, remove .cpp file
 */
class ProlongationPoint{

	unsigned int elm_1D_idx;
	unsigned int elm_3D_idx;
	unsigned int elm_3D_idx_old;

public:

	inline ProlongationPoint(){};
	inline ProlongationPoint(unsigned int elm1D, unsigned int elm3D, unsigned int elm3Dold):
			elm_1D_idx(elm1D), elm_3D_idx(elm3D), elm_3D_idx_old(elm3Dold){};
	inline ~ProlongationPoint(){};

	inline unsigned int get_elm_1D_idx() const{
		return elm_1D_idx;
	};

	inline unsigned int get_elm_3D_idx() const{
		return elm_3D_idx;
	};

	inline unsigned int get_elm_3D_idx_old() const{
		return elm_3D_idx_old;
	};
};

} // END NAMESPACE
#endif /* PROLONGATIONPOINT */
