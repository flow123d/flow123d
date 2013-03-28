/*
 * ProlongationPoint.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */
/*	SideIter elm_side = ele->side(1);
		Edge * edg = elm_side->edge();
		for(unsigned int i=0; i < edg->n_sides;i++) {
		SideIter other_side=edg->side(i);
		if (other_side != elm_side) {
		inspect.Add_element(idx);
		}
	}
	* Najde další 3D element, zkontrolu, jestli neni uz projity
	* Spočte průsečík s 1D PP, uloží průsečík, uloží průnik
	* Vytvoří nový 3D PP, eventuelně 1D PP, uloží jej do fronty 1D/3D a označí jej jako projité
	*/

class ProlongationPoint {
	 /// Intersection of prolongation point index
	    unsigned int intersect_idx;
	    /// Side of the intersection line
	    unsigned int local_side_1D;
	    /// Side of the 3D element of the intersection line (only for 3D prolongation points)
	    int local_side_3D;

	    /**
		Side stena = mesh->element[ intersect_.get_higher_dim_ele_idx() ].side( local_side_3D );
		**/
public:
	inline ProlongationPoint(unsigned int inter_idx, unsigned int local_side_1, int local_side_3):intersect_idx(inter_idx),local_side_1D(local_side_1),local_side_3D(local_side_3){};
	inline ~ProlongationPoint(){};
};
