/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "inspectelements.h"
namespace computeintersection {

InspectElements::InspectElements(){};

InspectElements::InspectElements(Mesh* _mesh):mesh(_mesh){};

InspectElements::~InspectElements(){};

void InspectElements::computeIntersections2d3dInit(){

	flag_for_3D_elements.assign(mesh->n_elements(), -1);

	FOR_ELEMENTS(mesh, elm) {
		if (elm->dim() == 2) {
			intersection_list[elm.index()] = std::vector<IntersectionLocal>();
			closed_elements[elm.index()] = false;
		}
	}

};

void InspectElements::computeIntersections2d3dUseProlongationTable(std::vector<std::pair<unsigned int, unsigned int>> &prolongation_table, const ElementFullIter &elm, const ElementFullIter &ele){
	//xprintf(Msg, "========PRODLUZUJI========\n");
	for(unsigned int i = 0; i < prolongation_table.size();i++){

		unsigned int stena = std::get<0>(prolongation_table[i]);
		unsigned int typ_elm = std::get<1>(prolongation_table[i]);

		if(typ_elm == 0){
			// prodlužuji hranou
			//xprintf(Msg,"Procházím hranu(%d) na id elementu(%d), hrana(%d)\n"
			//					, stena,elm->index(),elm->side(2-stena)->el_idx());

			SideIter elm_side = elm->side((unsigned int)(2-stena));
			Edge *edg = elm_side->edge();

			for(int j=0; j < edg->n_sides;j++) {
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					//xprintf(Msg, "\t\t Idx původního elementu a jeho hrany(%d,%d) - Idx sousedního elementu a jeho hrany(%d,%d)",elm->index(),stena,other_side->element()->index(),other_side->el_idx());
					unsigned int sousedni_element = other_side->element()->index(); // 2D element
						if(!intersectionExists(sousedni_element,ele->index())){
							//flag_for_3D_elements[ele->index()] = sousedni_element;

							// Vytvoření průniku bez potřeby počítání
							IntersectionLocal il_other(sousedni_element, ele->index());
							intersection_list[sousedni_element].push_back(il_other);

							ProlongationLine pl2(sousedni_element, ele->index(), intersection_list[sousedni_element].size() - 1);
							prolongation_line_queue_2D.push(pl2);
						}
				}
			}

		}else{
			// prodlužuji stěnou
			//xprintf(Msg,"Procházím stěnu(%d) na id elementu(%d), stěna(%d)\n"
			//		, stena,ele->index(),ele->side(3-stena)->el_idx());


			SideIter elm_side = ele->side((unsigned int)(3-stena)); //ele->side(3-stena);
			Edge *edg = elm_side->edge();

			for(int j=0; j < edg->n_sides;j++) {
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					//xprintf(Msg, "\t\t Idx původního elementu a jeho stěny(%d,%d) - Idx sousedního elementu a jeho stěny(%d,%d)\n",ele->index(),stena,other_side->element()->index(),other_side->el_idx());


					unsigned int sousedni_element = other_side->element()->index();
					if(flag_for_3D_elements[sousedni_element] == -1 || (flag_for_3D_elements[sousedni_element] != (int)elm->index() && !intersectionExists(elm->index(),sousedni_element))){
						flag_for_3D_elements[sousedni_element] = elm->index();
						// Jedná se o vnitřní čtyřstěny v trojúhelníku

						// Vytvoření průniku bez potřeby počítání
						IntersectionLocal il_other(elm.index(), sousedni_element);
						intersection_list[elm.index()].push_back(il_other);

						ProlongationLine pl2(elm.index(), sousedni_element, intersection_list[elm.index()].size() - 1);
						prolongation_line_queue_3D.push(pl2);

					}
				}
			}

		}

	}
};

void InspectElements::computeIntersections2d3dProlongation(const ProlongationLine &pl){

	// Měly bychom být vždy ve stejném trojúhelníku, tudíž updateTriangle by měl být zbytečný
	ElementFullIter elm = mesh->element(pl.getElement2DIdx());
	ElementFullIter ele = mesh->element(pl.getElement3DIdx());

	this->UpdateTriangle(elm);
	this->UpdateTetrahedron(ele);

	ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
	CI_23.init();
	CI_23.compute(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()]);

	std::vector<std::pair<unsigned int, unsigned int>> prolongation_table;
	intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].traceGenericPolygon(prolongation_table);

	if(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].getIPsize() > 2){
		all_intersections.push_back(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()]);
		computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);
	}
};


void InspectElements::ComputeIntersections23(){

	computeIntersections2d3dInit();

	unsigned int elementLimit = 20;
	BIHTree bt(mesh, elementLimit);

	FOR_ELEMENTS(mesh, elm) {
		if (elm->dim() == 2 && !closed_elements[elm.index()]) {

			this->UpdateTriangle(elm);


			std::vector<unsigned int> searchedElements;
			TTriangle tt;
			tt.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
						 TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
						 TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)) );

			//FieldInterpolatedP0<3,FieldValue<3>::Scalar>::create_triangle(efi,tt);
			BoundingBox elementBoundingBox = tt.get_bounding_box();
			bt.find_bounding_box(elementBoundingBox, searchedElements);
			bool prunik = false;

			for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
				int idx = *it;
				ElementFullIter ele = mesh->element( idx );
				if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] != (int)elm->index()) {

					this->UpdateTetrahedron(ele);

					IntersectionLocal il(elm.index(), ele->index());
					ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
					CI_23.init();
					CI_23.compute(il);
					//il.tracePolygon();
					std::vector<std::pair<unsigned int, unsigned int>> prolongation_table;
					il.traceGenericPolygon(prolongation_table);


					//xprintf(Msg, "Polygon(%d) - patological: %d \n",il.getIPsize(), il.isPatological());

					if(il.getIPsize() > 2){
						all_intersections.push_back(il);

						intersection_list[elm.index()].push_back(il);
						flag_for_3D_elements[ele->index()] = elm.index();
						xprintf(Msg,"Velikost intersection_list(%d), pocet IL v konkretnim listu(%d)\n", intersection_list.size(), intersection_list[elm.index()].size());

						// PRODLUŽOVÁNÍ
						computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);

						// Dokončila se trasovací tabulka = fronty jsou naplněny

						while(1){
						// Zpracování front
							while(!prolongation_line_queue_3D.empty()){

								computeIntersections2d3dProlongation(prolongation_line_queue_3D.front());
								prolongation_line_queue_3D.pop();

							}
							// Pridat priznak trojuhleniku, ze je projity
							if(!prolongation_line_queue_2D.empty()){
								computeIntersections2d3dProlongation(prolongation_line_queue_2D.front());
								prolongation_line_queue_2D.pop();
							}

							if(prolongation_line_queue_2D.empty() && prolongation_line_queue_3D.empty()){
								break;
							}

						}
						prunik = true;
						break; // ukončí procházení dalších bounding boxů
					}

				}
			}
			// Prošlo se celé pole sousedním bounding boxů, pokud nevznikl průnik, může se trojúhelník nastavit jako uzavřený
			closed_elements[elm.index()] = true;
			if(prunik){
				break; // do budoucna odstranit break - tady ho mam, abych měl jistotu, že se provádí průchod
			}
		}
	}

	/*xprintf(Msg,"Element(3616) - velikost(%d)\n", intersection_list[3616].size());
	for(unsigned int i = 0;i < intersection_list[3616].size();i++){
		IntersectionLocal il = intersection_list[3616][i];
		xprintf(Msg,"elementy(%d,%d) - obsah(%f), bodu(%d)\n", il.idx_2D(), il.idx_3D(), il.getArea(), il.getIPsize());
		il.printTracingTable();
	}
	xprintf(Msg,"Element(3617) - velikost(%d)\n", intersection_list[3617].size());*/
};


bool InspectElements::intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx){

	for(unsigned int i; i < intersection_list[elm_2D_idx].size();i++){

		if(intersection_list[elm_2D_idx][i].idx_3D() == elm_3D_idx){
			return true;
		}

	}
	return false;
};


void InspectElements::UpdateTriangle(const ElementFullIter &element_2D){

			arma::vec3 pole_bodu[3] = {element_2D->node[0]->point(),element_2D->node[1]->point(),element_2D->node[2]->point()};
		    triangle = Simplex<2>(pole_bodu);

};

void InspectElements::UpdateTetrahedron(const ElementFullIter &element_3D){

			arma::vec3 pole_bodu[4] = {element_3D->node[0]->point(),element_3D->node[1]->point(),element_3D->node[2]->point(),element_3D->node[3]->point()};
		    tetrahedron = Simplex<3>(pole_bodu);

};

void InspectElements::print(unsigned int vyber){


	FILE * soubor;
	if(vyber == 0){
		soubor = fopen("pokusTriangle.msh","w");

	}else{
		soubor = fopen("pokusTetrahedron.msh","w");

	}

	fprintf(soubor, "$MeshFormat\n");
	fprintf(soubor, "2.2 0 8\n");
	fprintf(soubor, "$EndMeshFormat\n");
	fprintf(soubor, "$Nodes\n");



	int pocet_pruseciku = 0;
	int pocet_predtim = 0;

	FOR_NODES(mesh, nod){
		pocet_predtim++;
	}

		for(unsigned int i = 0; i < all_intersections.size(); i++){
			pocet_pruseciku += all_intersections[i].getIPsize();
		}
	fprintf(soubor, "%d\n", (pocet_pruseciku + pocet_predtim));

	unsigned int velikost = 0;
	int pocet_prus = pocet_predtim;
	unsigned int velikost_pruseciku = all_intersections.size();

	FOR_NODES(mesh, nod){
		arma::vec3 bod = nod->point();
		fprintf(soubor,"%d %f %f %f\n", nod.id(), bod[0], bod[1], bod[2]);

	}

	for(unsigned int i = 0; i < velikost_pruseciku; i++){

		IntersectionLocal il = all_intersections[i];
		velikost = il.getIPsize();

			unsigned int idx2D = all_intersections[i].idx_2D();
			unsigned int idx3D = all_intersections[i].idx_3D();

			ElementFullIter el2D = mesh->element(idx2D);
			ElementFullIter el3D = mesh->element(idx3D);


			for(unsigned int j = 0; j < velikost; j++){
				pocet_prus++;

				IntersectionPoint<2,3> IP23 = il.get_point(j);
				arma::vec3 T_globalni;
				if(vyber == 0){
					T_globalni = (IP23.getLocalCoords1())[0] * el2D->node[0]->point()
													   +(IP23.getLocalCoords1())[1] * el2D->node[1]->point()
													   +(IP23.getLocalCoords1())[2] * el2D->node[2]->point();
				}else{
					T_globalni = (IP23.getLocalCoords2())[0] * el3D->node[0]->point()
													   +(IP23.getLocalCoords2())[1] * el3D->node[1]->point()
													   +(IP23.getLocalCoords2())[2] * el3D->node[2]->point()
													   +(IP23.getLocalCoords2())[3] * el3D->node[3]->point();
				}
				fprintf(soubor,"%d %f %f %f\n", pocet_prus, T_globalni[0], T_globalni[1], T_globalni[2]);
				//fprintf(soubor,"%d %f %f %f\n", pocet_prus, C_globalni[0], C_globalni[1], C_globalni[2]);

			}
	}




	fprintf(soubor,"$EndNodes\n");
	fprintf(soubor,"$Elements\n");
	fprintf(soubor,"%d\n", (pocet_pruseciku + mesh->n_elements()) );

	FOR_ELEMENTS(mesh, elee){
		// 1 4 2 30 26 1 2 3 4
		// 2 2 2 2 36 5 6 7
		if(elee->dim() == 3){
			int id1 = mesh->node_vector.index(elee->node[0]) + 1;
			int id2 = mesh->node_vector.index(elee->node[1]) + 1;
			int id3 = mesh->node_vector.index(elee->node[2]) + 1;
			int id4 = mesh->node_vector.index(elee->node[3]) + 1;

			fprintf(soubor,"%d 4 2 30 26 %d %d %d %d\n", elee.id(), id1, id2, id3, id4);
		}else if(elee->dim() == 2){
			int id1 = mesh->node_vector.index(elee->node[0]) + 1;
			int id2 = mesh->node_vector.index(elee->node[1]) + 1;
			int id3 = mesh->node_vector.index(elee->node[2]) + 1;
			fprintf(soubor,"%d 2 2 2 36 %d %d %d\n", elee.id(), id1, id2, id3);

		}else{

		}
		}


	pocet_prus = mesh->n_elements();
	int posledni = 0;
	int vrcholy = pocet_predtim;
		for(unsigned int i = 0; i < all_intersections.size(); i++){


			IntersectionLocal il = all_intersections[i];
			velikost = il.getIPsize();

				for(unsigned int j = 0; j < velikost; j++){
					pocet_prus++;
					vrcholy++;
					if(j == 0){
						posledni = vrcholy;
					}

					if((j+1) == velikost){
						fprintf(soubor,"%d 1 2 18 7 %d %d\n", pocet_prus, vrcholy, posledni);
					}else{
						fprintf(soubor,"%d 1 2 18 7 %d %d\n", pocet_prus, vrcholy, vrcholy+1);
					}

				}
		}

	fprintf(soubor,"$EndElements\n");
	fclose(soubor);

}

} // END namespace
