/*
 * ProlongationPoint3D.h
 *
 *  Created on: 20.3.2013
 *      Author: viktor
 */

class PPBase;
//class Intersection0D;

class ProlongationPoint3D : public PPBase{
	//Intersection0D *previous;
	//Side stena;
public:
	inline ProlongationPoint3D(){};
	inline ~ProlongationPoint3D(){};
	inline void ZjistiDalsi3DEl(){
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

	};
};


