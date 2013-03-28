/*
 * inspectelements.cpp
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 */

#include "fast_intersection/inspectelements.h"


InspectElements::InspectElements(){

}

InspectElements::InspectElements(Mesh* sit_){
	sit = sit_;
	projeti(sit->n_elements(),false);
	calculate_intersections();
}

void InspectElements::calculate_intersections(){
	unsigned int elementLimit = 20;
	GmshMeshReader reader(sit);
	reader.read_mesh(&sit);
	BIHTree bt(&sit, elementLimit);

	FOR_ELEMENTS(&sit, elm) {
		if (elm->dim() == 1 && !projeti[elm->index()]) {
			std::vector<unsigned int> searchedElements;
		    TAbscissa ta;
		    ElementFullIter efi = sit.element(elm.index());
		    FieldInterpolatedP0<3,FieldValue<3>::Scalar>::createAbscissa(efi, ta);
		    BoundingBox elementBoundingBox = ta.get_bounding_box();
		    bt.find_bounding_box(elementBoundingBox, searchedElements);

		    for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
		    	int idx = *it;
		        ElementFullIter ele = sit.element( idx );
		        if (ele->dim() == 3) {
		           	/* Spočte se průnik mezi 3D a 1D elementem, pokud jsou nalezeny průsečíky
		           	 * vytvoří se prolongation point a uloží se do fronty
		           	 * */


		        	Intersection_Local il; // vytvorit průnik
		        	all_intersection.push_back(il);

		        	// Pokud bude 1D protinat 3D cely -> vzniknou 2 Prolongation pointy
		        	// ! ošetřit všechny možné případy
		        	ProlongationPoint pp1(il.getID(),1, il.getSide(0));
		        	ProlongationPoint pp2(il.getID(),1, il.getSide(1));

		        	ppoint.push(pp1);
		        	ppoint.push(pp2);



		        	while(!ppoint.empty()){
		        		ppoint.front();

		        		// Výpočetní postup pro další průniky a elementy

		        		ppoint.pop();
		        	}
		        }
		    }
		}
	}
}

InspectElements::~InspectElements(){}

vector<Intersection_Local> InspectElements::getIntersections(){return all_intersection;}
