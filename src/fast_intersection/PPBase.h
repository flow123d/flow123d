/*
 * PPBase.h
 *
 *  Created on: 20.3.2013
 *      Author: viktor
 */

class Intersection0D;

class PPBase {
	//Element elm; //nebo index elementu

    /// Intersection of prolongation point index
    unsigned int intersect_idx;
    /// Side of the intersection line
    unsigned int local_side_1D;
    /// Side of the 3D element of the intersection line (only for 3D prolongation points)
    int local_side_3D;

    /**
	Intersection0D prusecik; // U 1D elementu to bude vlastně počáteční a koncový bod
	Side stena = mesh->element[ intersect_.get_higher_dim_ele_idx() ].side( local_side_3D );
	**/
public:
	inline PPBase(){};
	inline ~PPBase(){};
};

