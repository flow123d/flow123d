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
	plucker_product.assign(6, NULL);
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
			xprintf(Msg, "Nalezen 1D element \n");
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
		        	xprintf(Msg, "Nalezen 3D element \n");
		        	if(calculate_prolongation_point(elm, ele)){
		        		while(!ppoint.empty()){
		        			xprintf(Msg, "==============\n PROCHÁZENÍ FRONTY \n");
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

	std::vector<double> coords_3D;
	double theta;
	bool orientace;

	for(unsigned int i = 0; i < 4; i++){
		if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
			if(theta > 0 && theta < 1){
				ProlongationPoint pp(element_1D->index(), element_3D->index(), i, coords_3D, theta, orientace);
				ppoint.push(pp);
				xprintf(Msg, "SPOČETL SE PPOINT \n");
				xprintf(Msg, "ID1D: %d ID3D %d STENA %d --- THETA %f \n", element_1D->index(), element_3D->index(), i, theta);
				SideIter elm_side = element_3D->side(i);
				Edge *edg = elm_side->edge();
				ASSERT(edg, "Null edge \n");
				for(unsigned int j=0; j < edg->n_sides;j++) {
					SideIter other_side=edg->side(j);
					if (other_side != elm_side) {
						xprintf(Msg, "SPOČETL SE SOUSEDNÍ PPOINT \n");
						ProlongationPoint pp2(element_1D->index(), other_side->element()->index() , other_side->el_idx() , coords_3D, theta, !orientace);
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
	//update_abscissa(sit->element(point.idx_elm1D()), true);
	update_abscissa(sit->element(point.idx_elm1D()), point.getOrientation());

	int stena;
	std::vector<double> coords_3D;
	double theta;
	bool orientace;
	/* Cyklus přes všechny stěny
	 * Pokud je stěna jiná než jaká je z prolongation pointu -> spočte se průsečík
	 * Pokud byl nalezen, porovná se, zda leží na úsečce
	 * Pokud leží, vytvoří se intersection_local -> naplní se; najde se sousedící element a vytvoří se další PPoint
	 * Pokud neleží, spočte se interpolace souřadnic, vytvoří se Intersection_local, přes stěnu 1D elementu se hledají další průniky
	 */

	for(unsigned int i = 0; i < 4; i++){
		if(i != point.idx_side3D()){
			xprintf(Msg,"Prochází ostatních stěn: %d \n", i);
			if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
				projeti[point.idx_elm1D()] = true;
				xprintf(Msg, "SPOČETL SE PRŮSEČÍK OD PPOINTU! \n");
				xprintf(Msg, "STENA %d --- THETA %f \n", i, theta);
				if(theta < 1){
					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(),point.local_coords_1D());
					/* Když původní 1D element je v opačné orientaci, uložim i jinak orientovanou thétu
					 * abscissa bude vždy ve správně orientaci
					 * */
					if(!point.getOrientation()){theta = 1 - theta;}
					il.add_local_coord(coords_3D, theta);
					all_intersection.push_back(il);

					/*řešit číslování stěn!! */
					SideIter elm_side = sit->element(point.idx_elm3D())->side(i);
					Edge *edg = elm_side->edge();
					for(unsigned int j=0; j < edg->n_sides;j++) {
						SideIter other_side=edg->side(j);
						if (other_side != elm_side) {
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(point.idx_elm1D(), other_side->element()->index(),other_side->el_idx(), coords_3D, theta, point.getOrientation());
							ppoint.push(pp);
						}
					}
				}else{
					xprintf(Msg, "INTERPOLACE \n");
					std::vector<double> local_3D_coords = local_vector_interpolation(point.local_coords_3D_ref(),coords_3D,point.local_coords_1D_if(),theta, 1);

					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(), point.local_coords_1D());

					unsigned int side;
					if(point.getOrientation()){
						side = 1;
					}else{
						side = 0;
					}
					il.add_local_coord(local_3D_coords, side);
					all_intersection.push_back(il);

					SideIter elm_side = sit->element(point.idx_elm1D())->side(side);
					Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							if(!projeti[other_side->element().index()]){
							calculate_intersection_from_1D(other_side->element().index(), point.idx_elm3D(), local_3D_coords);
							}
						}
					}
				}
				break;
			}
		}
	}
}

void InspectElements::calculate_intersection_from_1D(unsigned int idx_1D, unsigned int idx_3D, std::vector<double> &interpolated_3D_coords){

	double theta;
	std::vector<double> coords_3D;
	bool orientace;
	bool nalezeni = false;
	unsigned int stena;

	update_abscissa(sit->element(idx_1D), true);
	projeti[idx_1D] = true;

	for(unsigned int i = 0; i < 4; i++){
		if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
			if((theta < 0 && orientace) || (theta > 0 && !orientace)){
				nalezeni = true;
				stena = i;
				break;
			}
		}
	}
	if(nalezeni){
			if(theta < 1 && theta > 0){
				IntersectionLocal il(idx_1D, idx_3D);
				//souřadnice z interpolaci, lokální théta = 0;
				unsigned int pocatecni_bod;
				if(orientace){
					pocatecni_bod = 1;
				}
				else{
					pocatecni_bod = 0;
				}
				il.add_local_coord(interpolated_3D_coords, pocatecni_bod);
				il.add_local_coord(coords_3D, theta);
				all_intersection.push_back(il);

				SideIter elm_side = sit->element(idx_3D)->side(stena);
				Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(idx_1D, other_side->element()->index(),other_side->el_idx(), coords_3D, theta, !orientace);
							ppoint.push(pp);
						}
					}

			}else{

				unsigned int pocatecni_bod;
				if(orientace){
					pocatecni_bod = 1;
				}
				else{
					pocatecni_bod = 0;
				}
				unsigned int koncovy_bod = (pocatecni_bod+1)%2;
				std::vector<double> local_3D_coords = local_vector_interpolation(interpolated_3D_coords,coords_3D,pocatecni_bod,theta, koncovy_bod);

				IntersectionLocal il(idx_1D, idx_3D);
				il.add_local_coord(interpolated_3D_coords, pocatecni_bod);
				il.add_local_coord(local_3D_coords, koncovy_bod);
				all_intersection.push_back(il);

									SideIter elm_side = sit->element(idx_1D)->side(koncovy_bod);
									Edge *edg = elm_side->edge();
									for(unsigned int j = 0; j < edg->n_sides; j++){
										SideIter other_side = edg->side(j);
										if(other_side != elm_side){
											if(!projeti[other_side->element().index()]){
											calculate_intersection_from_1D(other_side->element().index(), idx_3D, local_3D_coords);
											}
										}
									}

			}
	}


}

InspectElements::~InspectElements(){}

vector<IntersectionLocal> InspectElements::getIntersections(){return all_intersection;}

void InspectElements::update_tetrahedron(const ElementFullIter &element_3D){

			for(unsigned int i = 0; i < 6;i++){
				plucker_product[i] = NULL;
			}

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

void InspectElements::fill_plucker_product(int index1, int index2, int index3, double &c, double &d, double &e, Simplex<2,3> sm, int &stena){
	if(plucker_product[index1] == NULL){

		c = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint());
		if(stena == 2){
			c = -c;
		}
		plucker_product[index1] = new double(c);

	}
	else{
		c = *plucker_product[index1];
	}

	if(plucker_product[index2] == NULL){
		d = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
		plucker_product[index2] = new double(d);}
	else{
		d = *plucker_product[index2];
	}

	if(plucker_product[index3] == NULL){
		e = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());

		if(stena == 3){
			e = -e;
		}

		plucker_product[index3] = new double(e);}
	else{
		e = *plucker_product[index3];
	}
}


bool InspectElements::intersection_1D_2D(Simplex<2,3> sm, int stena, std::vector<double> &coords_3D, double &local_abscissa, bool &orientace){
	double c,d,e;
	xprintf(Msg, "FUNCE INTERSECTION 1D 2D jede \n");

	if(stena == 0){
		fill_plucker_product(0,1,2,c,d,e, sm, stena);
		xprintf(Msg, "C: %f D: %f E: %f \n", c,d,e);
	}
	else if(stena == 1){
		fill_plucker_product(0,3,4,c,d,e, sm, stena);


		//plucker_product[0] == NULL ? *plucker_product[0] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint());
		//plucker_product[3] == NULL ? *plucker_product[3] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint());
	    //plucker_product[4] == NULL ? *plucker_product[4] = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint());



		//c = *plucker_product[0]*(-1);
		//d = *plucker_product[3]*(-1);
		//e = *plucker_product[4]*(-1);
		c *= -1;
		d *= -1;
		e *= -1;
		xprintf(Msg, "C: %f D: %f E: %f \n", c,d,e);

	}
	else if(stena == 2){
		fill_plucker_product(2,5,4,c,d,e, sm, stena);
		//c = plucker_product[2] == NULL ? *plucker_product[2] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint()) : *plucker_product[2];
		//d = plucker_product[5] == NULL ? *plucker_product[5] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint()) : *plucker_product[5];
		//e = plucker_product[4] == NULL ? *plucker_product[4] = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint()) : *plucker_product[4];

		c *= -1;
		xprintf(Msg, "C: %f D: %f E: %f \n", c,d,e);
	}
	else{
		fill_plucker_product(1,5,3,c,d,e, sm, stena);
		//c = plucker_product[1] == NULL ? *plucker_product[1] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint()) : *plucker_product[1];
		//d = plucker_product[5] == NULL ? *plucker_product[5] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint()) : *plucker_product[5];
		//e = plucker_product[3] == NULL ? *plucker_product[3] = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint()) : *plucker_product[3];
		e *= -1;
		xprintf(Msg, "C: %f D: %f E: %f \n", c,d,e);
	}
	xprintf(Msg, "FUNCE INTERSECTION 1D 2D spočteny pluckery \n");

	/* Při výpočtu jsou P. souřadnice zorientovány, aby vycházely kladné hodnoty, pokud úsečka vstupuje do 4 stěnu
	 * a záporné pokud vystupuje
	 * */

	if((c > 0 && d > 0 && e > 0) || (c < 0 && d < 0 && e < 0)){
		// c = w0; d = w1; e = w2
		// lokální alfa = w2/soucet; lokální beta = w1/soucet; => lokální souřadnice na stěně
		double alfa = e/(c+d+e);
		double beta = c/(c+d+e);
		// lokální souřadnice na přímce T
		// T = localAbscissa= - A(i) + ( 1 - alfa - beta ) * V0(i) + alfa * V1(i) + beta * V2 (i) / U(i)
		// i = max z U(i)
		Vector<3> vec(abscissa.getPointA(),abscissa.getPointB());
		int i = 0;
		double max = vec[0];
		if(vec[1] > max) max = vec[1]; i = 1;
		if(vec[2] > max) max = vec[2]; i = 2;

		local_abscissa = (-abscissa.getPointA()[i] + (1 - alfa - beta)*sm[0][0].getPoint()[i] +
				alfa*sm[0][1].getPoint()[i] + beta*sm[1][1].getPoint()[i])/max;

		coords_3D.push_back(alfa); coords_3D.push_back(beta);

		if(c*d*e > 0){
			orientace = true; // orientace je v pořádku -> úsečka vchází do 4stěnu
		}
		else{
			orientace = false; // úsečka vychází
		}


		return true;
	}
	else{
		return false;
    }
};


} // namespace fast_1_3 close
