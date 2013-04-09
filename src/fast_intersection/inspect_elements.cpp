/*
 * inspectelements.cpp
 *
 *  Created on: 11.3.2013
 *      Author: viktor
 */

#include "system/system.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "fields/field_interpolated_p0.hh"
#include "system/sys_profiler.hh"

#include "inspect_elements.h"
//#include "bp_fris/intersections.h"

namespace fast_1_3{

InspectElements::InspectElements(){}

InspectElements::InspectElements(Mesh* sit_):sit(sit_){
	projeti.assign(sit->n_elements(),false);
	calculate_intersections();
}

double InspectElements::get_local_coords_1D(SPoint<3> a,SPoint<3> b,SPoint<3> x){
	Vector<3> vec(a,b);
	int i = 0;
	double max = vec[0];
	if(vec[1] > max) max = vec[1]; i = 1;
	if(vec[2] > max) max = vec[2]; i = 2;

	return ((x[i] - a[i]) / max);
}

SPoint<3> InspectElements::local_interpolation(SPoint<3> &point_A,SPoint<3> &point_B,const double &t_A,const double &t_B,const double &t){

	// 3D Point A, B, X: X = (tb - t)A + (t - ta)B / (tb - ta)

	double koef1 = t_B - t;
	double koef2 = t - t_A;
	double koef3 = t_B - t_A;

	return SPoint<3>((koef1*point_A[0] + koef2*point_B[0])/koef3,
						 (koef1*point_A[1] + koef2*point_B[1])/koef3,
						 (koef1*point_A[2] + koef2*point_B[2])/koef3);
}

std::vector<double> InspectElements::local_vector_interpolation(const std::vector<double> &point_A,const std::vector<double> &point_B,
		const double &t_A, const double &t_B, const double &t){

	double koef1 = t_B - t;
	double koef2 = t - t_A;
	double koef3 = t_B - t_A;
	std::vector<double> result;

	result.push_back((koef1*point_A[0] + koef2*point_B[0])/koef3);
	result.push_back((koef1*point_A[1] + koef2*point_B[1])/koef3);

	double temp = point_A.size() == 2 ? 1 - point_A[0] - point_A[1] : point_A[2];
	double temp2 = point_B.size() == 2 ? 1 - point_B[0] - point_B[1] : point_B[2];

	result.push_back((koef1*temp + koef2*temp2)/koef3);

	return result;
}

void InspectElements::calculate_intersections(){
	unsigned int elementLimit = 20;
	BIHTree bt(sit, elementLimit);

	FOR_ELEMENTS(sit, elm) {
		if (elm->dim() == 1 && !projeti[elm->index()]) {
			projeti[elm->index()] = true;
			std::vector<unsigned int> searchedElements;
		    TAbscissa ta;
		    ElementFullIter efi = sit->element(elm.index());
		    FieldInterpolatedP0<3,FieldValue<3>::Scalar>::createAbscissa(efi, ta);
		    BoundingBox elementBoundingBox = ta.get_bounding_box();
		    bt.find_bounding_box(elementBoundingBox, searchedElements);

		    for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
		    	int idx = *it;
		        ElementFullIter ele = sit->element( idx );
		        if (ele->dim() == 3) {
		        	if(calculate_prolongation_point(elm, ele)){
		        		while(!ppoint.empty()){
		        			calculate_from_prolongation_point(ppoint.front());
		        			ppoint.pop();
		        	    }
		        		break;
		        	}
		        }
		    }
		}
	}
}

bool InspectElements::calculate_prolongation_point(const ElementFullIter &element_1D, const ElementFullIter &element_3D){

	update_tetrahedron(element_3D);
	update_abscissa(element_1D, true);

	int stena = -1;
	std::vector<double> coords_3D;
	double theta;

	for(unsigned int i = 0; i < 4; i++){
		if(IntersectionsOp_1D_2D(abscissa, tetrahedron[i], i, coords_3D, theta)){
			if(theta > 0 && theta < 1){
				ProlongationPoint pp(element_1D->index(), element_3D->index(), i, coords_3D, theta);
				ppoint.push(pp);

				SideIter elm_side = element_3D->side(i);
				Edge *edg = elm_side->edge();
				ASSERT(edg, "Null edge \n");
				for(unsigned int j=0; j < edg->n_sides;j++) {
					SideIter other_side=edg->side(j);
					if (other_side != elm_side) {
						// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
						// možná řešit i už spočítané pluckerovy souřadnice
						stena = i; // doplnit
						//sit->element(other_side->element())
						ProlongationPoint pp2(element_1D->index(), other_side->element()->index() , other_side->el_idx() , coords_3D, theta);
						ppoint.push(pp2);
					}
				}
				return true;
			}
		}
	}
	return false;
}

void InspectElements::calculate_from_prolongation_point(ProlongationPoint &point){

	update_tetrahedron(sit->element(point.idx_elm3D()));
	// orientace z PPointu
	update_abscissa(sit->element(point.idx_elm1D()), true);

	//zjistit číslování stěn, správně permutovat
	int stena;
	std::vector<double> coords_3D;
	double theta;
	/* Cyklus přes všechny stěny
	 * Pokud je stěna jiná než jaká je z prolongation pointu -> spočte se průsečík
	 * Pokud byl nalezen, porovná se, zda leží na úsečce
	 * Pokud leží, vytvoří se intersection_local -> naplní se; najde se sousedící element a vytvoří se další PPoint
	 * Pokud neleží, spočte se interpolace souřadnic, vytvoří se Intersection_local, přes stěnu 1D elementu se hledají další průniky
	 */

	for(unsigned int i = 0; i < 4; i++){
		if(i != point.idx_side3D()){
			if(IntersectionsOp_1D_2D(abscissa, tetrahedron[i], i, coords_3D, theta)){
				// zajistit, aby úsečka byla ve směru 0 -> 1
				projeti[point.idx_elm1D()] = true;
				if(theta < 1){
					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(),point.local_coords_1D());
					il.add_local_coord(coords_3D, theta);
					all_intersection.push_back(il);

					/*řešit číslování stěn!! */
					SideIter elm_side = sit->element(point.idx_elm3D())->side(i);
					Edge *edg = elm_side->edge();
					for(unsigned int j=0; j < edg->n_sides;j++) {
						SideIter other_side=edg->side(j);
						if (other_side != elm_side) {
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(point.idx_elm1D(), other_side->element()->index(),other_side->el_idx(), coords_3D, theta);
							ppoint.push(pp);
						}
					}
				}else{
					// ověřit zda je orientace 1D elementu ve směru 0 -> 1
					// Interpolation:
					std::vector<double> local_3D_coords = local_vector_interpolation(point.local_coords_3D_ref(),coords_3D,point.local_coords_1D(),theta, 1);

					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(), point.local_coords_1D());
					il.add_local_coord(local_3D_coords, 1);
					all_intersection.push_back(il);

					SideIter elm_side = sit->element(point.idx_elm1D())->side(1);
					Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							if(!projeti[other_side->element().index()]){
							calculate_intersection_from_1D(other_side->element().index(), point.idx_elm3D());
							}
						}
					}
				}
				break;
			}
		}
	}
}

void InspectElements::calculate_intersection_from_1D(unsigned int idx_1D, unsigned int idx_3D){

	double theta;
	std::vector<double> coords_3D;

	// orientace! koncový bod úsečky, je počáteční nové
	update_abscissa(sit->element(idx_1D), true);
	projeti[idx_1D] = true;

	// Opět ověřit orientaci úsečky ve směru 0->1
	for(unsigned int i = 0; i < 4; i++){
		if(IntersectionsOp_1D_2D(abscissa, tetrahedron[i], i, coords_3D, theta)){
			// přidat break; zpracování až za for cyklem
			if(theta < 1 && theta > 0){
				IntersectionLocal il(idx_1D, idx_3D);
				//souřadnice z interpolaci, lokální théta = 0;

				// souřadnice z průsečíku
				il.add_local_coord(coords_3D, theta);
				all_intersection.push_back(il);

				SideIter elm_side = sit->element(idx_3D)->side(i);
				Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(idx_1D, other_side->element()->index(),other_side->el_idx(), coords_3D, theta);
							ppoint.push(pp);
						}
					}
				break;
			}else if(theta > 1){
				// úsečka končí ve vnitř -> procházet přes stěny všechny neprojité úsečky a metodu opakovat

				break;
			}


		}
	}

}

InspectElements::~InspectElements(){}

vector<IntersectionLocal> InspectElements::getIntersections(){return all_intersection;}

void InspectElements::update_tetrahedron(const ElementFullIter &element_3D){

			SPoint<3> spsim1; spsim1.setCoords(element_3D->node[0]->point());
			SPoint<3> spsim2; spsim2.setCoords(element_3D->node[1]->point());
			SPoint<3> spsim3; spsim3.setCoords(element_3D->node[2]->point());
			SPoint<3> spsim4; spsim4.setCoords(element_3D->node[3]->point());
		    SPoint<3> pole_bodu[4] = {spsim1,spsim2,spsim3,spsim4};
			tetrahedron = Simplex<3,3>(pole_bodu);
};

void InspectElements::update_abscissa(const ElementFullIter &element_1D, bool orientace){
	SPoint<3> sphp1;
	SPoint<3> sphp2;

	if(orientace){
    sphp1.setCoords(element_1D->node[0]->point());
	sphp2.setCoords(element_1D->node[1]->point());
	}else{
		sphp1.setCoords(element_1D->node[1]->point());
		sphp2.setCoords(element_1D->node[0]->point());
	}
	abscissa = HyperPlane<1,3>(sphp1, sphp2);
};


} // namespace fast_1_3 close
