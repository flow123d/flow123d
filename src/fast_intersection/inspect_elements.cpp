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

	Profiler::initialize(MPI_COMM_WORLD);
	{ START_TIMER("Test_bez_BIHtree");

	FOR_ELEMENTS(sit, elm) {
		if (elm->dim() == 1 && !projeti[elm->index()]) {
			xprintf(Msg, "-----Nalezen 1D element------ \n");
			//projeti[elm->index()] = true;
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
		        	//xprintf(Msg, "Nalezen 3D element \n");
		        	if(calculate_prolongation_point(elm, ele)){
		        		//xprintf(Msg, "Nalezen průsečík s 3D element \n");
		        		while(!ppoint.empty()){
		        			//xprintf(Msg, "PROCHÁZENÍ FRONTY %d\n", ppoint.size());
		        			calculate_from_prolongation_point(ppoint.front());
		        			ppoint.pop();
		        			//xprintf(Msg, "PPoint projety, velikost fronty: %d\n", ppoint.size());

		        	    }
		        		break;
		        	}
		        }
		    }
		}
	}
	END_TIMER("Test_bez_BIHtree");}
		Profiler::instance()->output(cout);
		Profiler::uninitialize();
}

bool InspectElements::calculate_prolongation_point(const ElementFullIter &element_1D, const ElementFullIter &element_3D){

	update_tetrahedron(element_3D);
	update_abscissa(element_1D, true);

	/*	Stěny elementu jsou jinak číslované než stěny mých simplexů!!! Všudě potřeba změnit
	 *  Stěna Simplexu | Stěna elementu
	 *   	0					3
	 *   	1					2
	 *   	2					1
	 *   	3					0
	 *
	 * */

	std::vector<double> coords_3D;
	double theta;
	bool orientace;

	for(unsigned int i = 0; i < 4; i++){
		if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
			if(theta >= 0 && theta <= 1){
				ProlongationPoint pp(element_1D->index(), element_3D->index(), i, coords_3D, theta, orientace);
				ppoint.push(pp);
				//xprintf(Msg, "CPP for THETA: %f\n", theta);
				//xprintf(Msg, "SPOČETL SE PPOINT \n");
				//xprintf(Msg, "ID1D: %d ID3D %d STENA %d --- THETA %f \n", element_1D->index(), element_3D->index(), i, theta);
				SideIter elm_side = element_3D->side(3-i); // Záměna stěn
				Edge *edg = elm_side->edge();
				//ASSERT(edg, "Null edge \n");
				for(unsigned int j=0; j < edg->n_sides;j++) {
					SideIter other_side=edg->side(j);
					if (other_side != elm_side) {
						//xprintf(Msg, "SPOČETL SE SOUSEDNÍ PPOINT \n");
						ProlongationPoint pp2(element_1D->index(), other_side->element()->index() , 3-other_side->el_idx(),
								coords_3D, theta, !orientace); // Záměna stěn
						//xprintf(Msg, "CPP for soused THETA: %f\n", theta);
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
	//xprintf(Msg,"CALCULATE_FROM_PROLONGATION_POINT:\n");
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
			//xprintf(Msg,"Prochází ostatních stěn: %d \n", i);
			if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
				projeti[point.idx_elm1D()] = true;

				//xprintf(Msg, "STENA %d --- THETA %f \n", i, theta);
				if(theta <= 1 && theta >= 0){
					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(),point.local_coords_1D());
					/* Když původní 1D element je v opačné orientaci, uložim i jinak orientovanou thétu
					 * abscissa bude vždy ve správně orientaci
					 * */
					if(!point.getOrientation()){theta = 1 - theta;}
					il.add_local_coord(coords_3D, theta);
					//xprintf(Msg, "CFPP for THETA: %f\n", theta);
					all_intersection.push_back(il);
					xprintf(Msg, "PRŮNIK! počet: %d\n", all_intersection.size());
					//xprintf(Msg, "VYTVOŘEN PRŮNIK(odppointu): THETA ALFA: %f THETA BETA: %f \n", point.local_coords_1D(), theta);

					/*řešit číslování stěn!! */
					SideIter elm_side = sit->element(point.idx_elm3D())->side(3-i); // Záměna stěn
					Edge *edg = elm_side->edge();
					for(unsigned int j=0; j < edg->n_sides;j++) {
						SideIter other_side=edg->side(j);
						if (other_side != elm_side) {
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(point.idx_elm1D(), other_side->element()->index(),
									3-other_side->el_idx(), coords_3D, theta, point.getOrientation()); // Záměna stěn
							//xprintf(Msg, "CFPP for soused THETA: %f\n", theta);
							ppoint.push(pp);
						}
					}
				}else{

					if(theta < 0){xprintf(Msg, "PROBLEM");
					point.getOrientation() ? xprintf(Msg, " orientace true\n") : xprintf(Msg, " orientace false\n");
					}
					//xprintf(Msg, "CFPP neni mezi 0 a 1 THETA: %f\n", theta);
					//xprintf(Msg, "INTERPOLACE \n");
					std::vector<double> local_3D_coords = local_vector_interpolation(point.local_coords_3D_ref(),coords_3D,point.local_coords_1D_if(),theta, 1);

					IntersectionLocal il(point.idx_elm1D(), point.idx_elm3D());
					il.add_local_coord(point.local_coords_3D(), point.local_coords_1D());

					double side;
					if(point.getOrientation()){
						side = 1.0;
					}else{
						side = 0.0;
					}
					il.add_local_coord(local_3D_coords, side);
					all_intersection.push_back(il);
					xprintf(Msg, "PRŮNIK! počet: %d\n", all_intersection.size());
					//xprintf(Msg, "VYTVOŘEN PRŮNIK(interpolace): THETA ALFA: %f THETA BETA: %f \n", point.local_coords_1D(), side);
					SideIter elm_side = sit->element(point.idx_elm1D())->side((unsigned int)side);
					Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							//xprintf(Msg, "Stěna %d \n", j);
							//if(true){
							if(!projeti[other_side->element().index()]){
								xprintf(Msg, "NALEZEN NEPROJITY 1D ELEMENT \n");
							calculate_intersection_from_1D(other_side->element().index(), point.idx_elm3D(), local_3D_coords, point.idx_elm1D(), point.getOrientation());
							}
						}
					}
				}
				break;
			}
		}
	}
}

void InspectElements::calculate_intersection_from_1D(unsigned int idx_1D, unsigned int idx_3D, std::vector<double> &interpolated_3D_coords, unsigned int idx_1D_previous, bool orientace_previous){

	double theta;
	std::vector<double> coords_3D;
	bool orientace;
	bool nalezeni = false;
	bool otoceni = false;
	unsigned int stena;

	ElementFullIter previous = sit->element(idx_1D_previous);
	ElementFullIter current = sit->element(idx_1D);
	if(orientace_previous){
		if(previous->node[1] == current->node[0]){
			otoceni = true;
		}
	}else{
		if(previous->node[0] == current->node[0]){
					otoceni = true;
		}
	}




	update_abscissa(sit->element(idx_1D), otoceni);


	for(unsigned int i = 0; i < 4; i++){

		if(intersection_1D_2D(tetrahedron[i], i, coords_3D, theta, orientace)){
			if(!orientace && theta >= 0){
				nalezeni = true;
				stena = i;
				projeti[idx_1D] = true;
				//xprintf(Msg, "Nalezeno v 1D\n");
				break;
			}
		}

	}
	if(nalezeni){
			if(theta <= 1 && theta >= 0){
				IntersectionLocal il(idx_1D, idx_3D);
				//souřadnice z interpolaci, lokální théta = 0;
				double pocatecni_bod;
				if(!otoceni){
					pocatecni_bod = 1.0;
				}
				else{
					pocatecni_bod = 0.0;
				}
				if(!otoceni){
					theta = 1 - theta;
				}
				il.add_local_coord(interpolated_3D_coords, pocatecni_bod);
				il.add_local_coord(coords_3D, theta);
				//xprintf(Msg, "F1D nalezeni for THETA: %f\n", theta);
				all_intersection.push_back(il);
				xprintf(Msg, "PRŮNIK! počet: %d\n", all_intersection.size());
				//xprintf(Msg, "VYTVOŘEN PRŮNIK(1D metoda): THETA ALFA: %f THETA BETA: %f \n", pocatecni_bod, theta);
				SideIter elm_side = sit->element(idx_3D)->side(3-stena); // Záměna stěn
				Edge *edg = elm_side->edge();
					for(unsigned int j = 0; j < edg->n_sides; j++){
						SideIter other_side = edg->side(j);
						if(other_side != elm_side){
							//xprintf(Msg, "vytvořen ppoint sousedící od 1 metody\n");
							//xprintf(Msg,"%d %d %d %f \n",idx_1D, other_side->element()->index(),other_side->el_idx(),theta);
							// řešit permutaci hran a jak jsou označené stěny!!! změnit alfu a betu
							ProlongationPoint pp(idx_1D, other_side->element()->index(),
									(3-other_side->el_idx()), coords_3D, theta, otoceni); // Záměna stěn
							ppoint.push(pp);
							//xprintf(Msg, "velikost fronty: %d \n", ppoint.size());
						}
					}

			}else{

				double pocatecni_bod;
				double koncovy_bod;
				if(!otoceni){
					pocatecni_bod = 1.0;
					koncovy_bod = 0;
				}
				else{
					pocatecni_bod = 0;
					koncovy_bod = 1.0;
				}

				std::vector<double> local_3D_coords = local_vector_interpolation(interpolated_3D_coords,coords_3D,pocatecni_bod,theta, koncovy_bod);

				IntersectionLocal il(idx_1D, idx_3D);
				il.add_local_coord(interpolated_3D_coords, pocatecni_bod);
				il.add_local_coord(local_3D_coords, koncovy_bod);
				all_intersection.push_back(il);
				xprintf(Msg, "PRŮNIK! počet: %d\n", all_intersection.size());
				//xprintf(Msg, "F1D intersection for THETA: %f -- ale pocatecni %f a koncovy %f\n", theta, pocatecni_bod, koncovy_bod);
				//xprintf(Msg, "VYTVOŘEN PRŮNIK(cela ve vnitr): THETA ALFA: %f THETA BETA: %f \n", pocatecni_bod, koncovy_bod);
									SideIter elm_side = sit->element(idx_1D)->side((unsigned int)koncovy_bod);
									Edge *edg = elm_side->edge();
									for(unsigned int j = 0; j < edg->n_sides; j++){
										//xprintf(Msg, "kolik edges: %d\n", edg->n_sides);
										SideIter other_side = edg->side(j);
										if(other_side != elm_side){
											//xprintf(Msg, "nalezena jina side\n");
											if(true){
											//if(!projeti[other_side->element().index()]){
											calculate_intersection_from_1D(other_side->element().index(), idx_3D, local_3D_coords, idx_1D, otoceni);
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

	for(unsigned int i = 0; i < 6;i++){
		plucker_product[i] = NULL;
	}

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

void InspectElements::fill_plucker_product(int index1, int index2, int index3, double &c, double &d, double &e, Simplex<2,3> &sm, int &stena){
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


bool InspectElements::intersection_1D_2D(Simplex<2,3> &sm, int stena, std::vector<double> &coords_3D, double &local_abscissa, bool &orientace){
	double c,d,e;

	if(stena == 0){
		fill_plucker_product(0,1,2,c,d,e, sm, stena);
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
	}
	else if(stena == 2){
		fill_plucker_product(2,5,4,c,d,e, sm, stena);
		//c = plucker_product[2] == NULL ? *plucker_product[2] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint()) : *plucker_product[2];
		//d = plucker_product[5] == NULL ? *plucker_product[5] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint()) : *plucker_product[5];
		//e = plucker_product[4] == NULL ? *plucker_product[4] = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint()) : *plucker_product[4];

		c *= -1;
	}
	else{
		fill_plucker_product(1,5,3,c,d,e, sm, stena);
		//c = plucker_product[1] == NULL ? *plucker_product[1] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint()), Vector<3>(sm[0][0].getPoint(),sm[0][1].getPoint())*(Vector<3>)sm[0][0].getPoint()) : *plucker_product[1];
		//d = plucker_product[5] == NULL ? *plucker_product[5] = abscissa.getPlucker()*Plucker(Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint()), Vector<3>(sm[0][1].getPoint(),sm[1][1].getPoint())*(Vector<3>)sm[0][1].getPoint()) : *plucker_product[5];
		//e = plucker_product[3] == NULL ? *plucker_product[3] = abscissa.getPlucker()*Plucker(Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint()), Vector<3>(sm[1][1].getPoint(),sm[0][0].getPoint())*(Vector<3>)sm[1][1].getPoint()) : *plucker_product[3];

		c *= -1;
		d *= -1;
		//e *= -1;
	}

	/* Při výpočtu jsou P. souřadnice zorientovány, aby vycházely kladné hodnoty, pokud úsečka vstupuje do 4 stěnu
	 * a záporné pokud vystupuje
	 * */
	//xprintf(Msg, "C: %f D: %f E: %f \n", c,d,e);

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
			orientace = false; // úsečka vychází
		}
		else{
			orientace = true; // orientace je v pořádku -> úsečka vchází do 4stěnu
		}



		return true;
	}
	else{

		return false;
    }
};

void InspectElements::print(char *nazev){

	FILE * soubor;
	soubor = fopen(nazev,"w");


	fprintf(soubor, "$MeshFormat\n");
	fprintf(soubor, "2.2 0 8\n");
	fprintf(soubor, "$EndMeshFormat\n");
	fprintf(soubor, "$Nodes\n");


	int aa = 2*all_intersection.size();

	//xprintf(Msg, "\n%d\n",aa);
	fprintf(soubor, "%d\n", aa);

	ElementFullIter ele1D = sit->element(0);

	for(unsigned int i = 0; i < all_intersection.size(); i++){

		IntersectionLocal il = all_intersection[i];
		ElementFullIter el1D = sit->element(il.idx_1D());

		double x = el1D->node[0]->getX() + (el1D->node[1]->getX() - el1D->node[0]->getX())*il.get_point(0)->el2_coord();
		double y = el1D->node[0]->getY() + (el1D->node[1]->getY() - el1D->node[0]->getY())*il.get_point(0)->el2_coord();
		double z = el1D->node[0]->getZ() + (el1D->node[1]->getZ() - el1D->node[0]->getZ())*il.get_point(0)->el2_coord();

		double x2 = el1D->node[0]->getX() + (el1D->node[1]->getX() - el1D->node[0]->getX())*il.get_point(1)->el2_coord();
		double y2 = el1D->node[0]->getY() + (el1D->node[1]->getY() - el1D->node[0]->getY())*il.get_point(1)->el2_coord();
		double z2 = el1D->node[0]->getZ() + (el1D->node[1]->getZ() - el1D->node[0]->getZ())*il.get_point(1)->el2_coord();

		fprintf(soubor,"%d %f %f %f\n", 2*i + 1, x, y, z);
		fprintf(soubor,"%d %f %f %f\n", 2*i + 2, x2, y2, z2);
		//fprintf(out_file,"%d %f %f %f\n", 2*i, x, y, z);
		//fprintf(out_file,"%d %f %f %f\n", 2*i + 1, x2, y2, z2);
	}

	int aaa = all_intersection.size();

	fprintf(soubor,"$EndNodes\n");
	fprintf(soubor,"$Elements\n");
	fprintf(soubor,"%d\n", aaa);

	for(int j = 0; j < aaa; j++){

		fprintf(soubor,"%d 1 2 18 7 %d %d\n", j+1, 2*j +1, 2*j +2);


	}
	fprintf(soubor,"$EndElements\n");
	fclose(soubor);

}



} // namespace fast_1_3 close
