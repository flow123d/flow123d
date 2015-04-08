/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "inspectelements.h"
namespace computeintersection {

InspectElements::InspectElements(Mesh* _mesh):mesh(_mesh){
	element_2D_index = -1;
};

InspectElements::~InspectElements(){};

template<>
void InspectElements::compute_intersections_init<1,3>(){
	flag_for_3D_elements.assign(mesh->n_elements(), -1);
	closed_elements.assign(mesh->n_elements(), false);
	intersection_line_list.assign(mesh->n_elements(),std::vector<IntersectionLine>());

	if(elements_bb.size() == 0){
		elements_bb.reserve(mesh->n_elements());
		bool first_3d_element = true;
		FOR_ELEMENTS(mesh, elm) {

			elements_bb[elm->index()] = elm->bounding_box();

				if (elm->dim() == 3){
					if(first_3d_element){
						first_3d_element = false;
						mesh_3D_bb = elements_bb[elm->index()];
					}else{
						mesh_3D_bb.expand(elements_bb[elm->index()].min());
						mesh_3D_bb.expand(elements_bb[elm->index()].max());
					}
				}
		}
	}
};

template<>
void InspectElements::compute_intersections_init<2,3>(){
	flag_for_3D_elements.assign(mesh->n_elements(), -1);
	closed_elements.assign(mesh->n_elements(), false);
	intersection_list.assign(mesh->n_elements(),std::vector<IntersectionLocal>());

	if(elements_bb.size() == 0){
		elements_bb.reserve(mesh->n_elements());
		bool first_3d_element = true;
		FOR_ELEMENTS(mesh, elm) {

			elements_bb[elm->index()] = elm->bounding_box();

				if (elm->dim() == 3){
					if(first_3d_element){
						first_3d_element = false;
						mesh_3D_bb = elements_bb[elm->index()];
					}else{
						mesh_3D_bb.expand(elements_bb[elm->index()].min());
						mesh_3D_bb.expand(elements_bb[elm->index()].max());
					}
				}
		}
	}
};

template<>
void InspectElements::compute_intersections<1,3>(){

	compute_intersections_init<1,3>();
	//computeIntersections2d3dInit();
	BIHTree bt(mesh, 20);

	FOR_ELEMENTS(mesh, elm) {
		if (elm->dim() == 1 && !closed_elements[elm->index()]&& elements_bb[elm->index()].intersect(mesh_3D_bb)) {
			xprintf(Msg, "-----Nalezen 1D element------ \n");

			update_abscissa(elm);
			std::vector<unsigned int> searchedElements;
			bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

			for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
				int idx = *it;
				ElementFullIter ele = mesh->element( idx );
				if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] == -1) {

					update_tetrahedron(ele);

					IntersectionLine il(elm->index(), ele->index());
					ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
					CI_13.init();
					CI_13.compute(il.get_points());

					if(il.size() > 1){
						closed_elements[elm->index()] = true;
						flag_for_3D_elements[ele->index()] = elm->index();
						intersection_line_list[elm->index()].push_back(il);

						prolongate_elements(il, elm, ele);

						while(!prolongation_point_queue.empty()){

							prolongate((ProlongationPoint)prolongation_point_queue.front());
							prolongation_point_queue.pop();

						}
						break;
					}
				}
			}
		}
	}
}

void InspectElements::prolongate_elements(const IntersectionLine &il, const ElementFullIter &elm, const ElementFullIter &ele){
	for(unsigned int i = 0; i < il.size();i++){

		if(il[i].isVertex()){

			SideIter elm_side = elm->side((unsigned int)(1-il[i].get_local_coords1()[0])); //ele->side(3-stena);
			Edge *edg = elm_side->edge();
			for(int j=0; j < edg->n_sides;j++) {

				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {

					unsigned int sousedni_element = other_side->element()->index();
					if(!closed_elements[sousedni_element]){
						prolongate_1D_element(other_side->element(), ele);
					}
				}
			}
		}else{

			SideIter elm_side = ele->side((unsigned int)(3-il[i].getSide2())); //ele->side(3-stena);
			Edge *edg = elm_side->edge();
			for(int j=0; j < edg->n_sides;j++) {

				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					unsigned int sousedni_element = other_side->element()->index();

					if(flag_for_3D_elements[sousedni_element] == -1){ // || (flag_for_3D_elements[sousedni_element] != (int)elm->index() && !intersectionExists(elm->index(),sousedni_element))){
						flag_for_3D_elements[sousedni_element] = elm->index();
						ProlongationPoint pp(elm->index(), sousedni_element, ele->index());
						prolongation_point_queue.push(pp);
					}
				}
			}
		}
	}
};

template<>
void InspectElements::compute_intersections<2,3>(){
	{ START_TIMER("Incializace pruniku");
		compute_intersections_init<2,3>();
		END_TIMER("Inicializace pruniku");}

		/*UpdateTriangle(mesh->element(5775));
		UpdateTetrahedron(mesh->element(172349));
		IntersectionLocal il(5775, 172349);
		ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
		CI_23.init();
		CI_23.compute(il);
		std::vector<unsigned int> prolongation_table;
		il.traceGenericPolygon(prolongation_table);
		intersection_list[5775].push_back(il);
		return;*/

		BIHTree bt(mesh, 20);
		bool nalezen = false;

		{ START_TIMER("Prochazeni vsech elementu");
		FOR_ELEMENTS(mesh, elm) {
			if (elm->dim() == 2 && !closed_elements[elm.index()] && elements_bb[elm->index()].intersect(mesh_3D_bb)) {
				update_triangle(elm);
				std::vector<unsigned int> searchedElements;
				bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

				{ START_TIMER("Hlavni vypocet");
				for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
				{
					unsigned int idx = *it;
					ElementFullIter ele = mesh->element(idx);

					if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] != (int)(elm->index())) {
						update_tetrahedron(ele);
						IntersectionLocal il(elm.index(), ele->index());
						ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
						CI_23.init();
						CI_23.compute(il);

						if(il.size() > 2){

							nalezen = true;
							std::vector<unsigned int> prolongation_table;
							//std::vector<unsigned int> pp_table;
							//il.trace_polygon_opt(pp_table);
							il.traceGenericPolygon(prolongation_table);
							//il.printTracingTable();

							{ START_TIMER("Prochazeni vsech front");

							intersection_list[elm.index()].push_back(il);
							flag_for_3D_elements[ele->index()] = elm.index();
							xprintf(Msg,"Velikost intersection_list(%d), pocet IL v konkretnim listu(%d)\n", intersection_list.size(), intersection_list[elm.index()].size());

							// PRODLUŽOVÁNÍ
							computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);

							while(1){
								//xprintf(Msg, "Fronta - prochazim(%d,%d)\n", elm->index(), ele->index());
							// Zpracování front
								while(!prolongation_line_queue_3D.empty()){
									//xprintf(Msg, "predwhile1\n");
									computeIntersections2d3dProlongation(prolongation_line_queue_3D.front());
									prolongation_line_queue_3D.pop();
									//xprintf(Msg, "while1\n");
								}

								if(element_2D_index >= 0){
									closed_elements[element_2D_index] = true;
								}

								element_2D_index = -1;
								// Pridat priznak trojuhleniku, ze je projity
								if(!prolongation_line_queue_2D.empty()){
									//xprintf(Msg, "predif1\n");
									computeIntersections2d3dProlongation(prolongation_line_queue_2D.front());
									prolongation_line_queue_2D.pop();
									//xprintf(Msg, "if1\n");
								}

								if(element_2D_index >= 0){
									closed_elements[element_2D_index] = true;
									//xprintf(Msg, "if2\n");
								}

								element_2D_index = -1;

								if(prolongation_line_queue_2D.empty() && prolongation_line_queue_3D.empty()){
									//xprintf(Msg, "if3\n");
									break;
								}

							}
							END_TIMER("Prochazeni vsech front");}
							//prunik = true;
							break; // ukončí procházení dalších bounding boxů
						}

					}
				}
				// Prošlo se celé pole sousedním bounding boxů, pokud nevznikl průnik, může se trojúhelník nastavit jako uzavřený
				closed_elements[elm.index()] = true;
				/*if(nalezen){
					break;
				}*/
				END_TIMER("Hlavni vypocet");}

			}
		}

		END_TIMER("Prochazeni vsech elementu");}

};

void InspectElements::prolongate_1D_element(const ElementFullIter &elm, const ElementFullIter &ele){

	update_abscissa(elm);
	closed_elements[elm->index()] = true;

	IntersectionLine il(elm->index(), ele->index());
	ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
	CI_13.init();
	CI_13.compute(il.get_points());

	if(il.size() > 1){
		intersection_line_list[elm->index()].push_back(il);
		prolongate_elements(il, elm, ele);
	}
};

void InspectElements::prolongate(const ProlongationPoint &pp){

	ElementFullIter elm = mesh->element(pp.get_elm_1D_idx());
	ElementFullIter ele = mesh->element(pp.get_elm_3D_idx());

	update_abscissa(elm);
	update_tetrahedron(ele);

	IntersectionLine il(elm->index(), ele->index());
	ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
	CI_13.init();
	CI_13.compute(il.get_points());

	if(il.size() > 1){
		intersection_line_list[elm->index()].push_back(il);
		prolongate_elements(il, elm, ele);
	}

}

void InspectElements::computeIntersections2d3dUseProlongationTable(std::vector<unsigned int> &prolongation_table, const ElementFullIter &elm, const ElementFullIter &ele){


	unsigned int elm_2D = ele->index();
	//xprintf(Msg, "========PRODLUZUJI========\n");
	for(unsigned int i = 0; i < prolongation_table.size();i++){

		unsigned int stena;
		unsigned int typ_elm;

		if(prolongation_table[i] >= 4){
			stena = prolongation_table[i] - 4;
			typ_elm = 0;
		}else{
			stena = prolongation_table[i];
			typ_elm = 1;
		}

		//cout << "[" << stena << "," << typ_elm << "]" << endl;


		if(typ_elm == 0){
			// prodlužuji hranou

			//xprintf(Msg,"Procházím hranu(%d) na id elementu(%d), hrana(%d)\n"
			//					, stena,elm->index(),elm->side(2-stena)->el_idx());


			SideIter elm_side = elm->side((3-stena)%3);


			Edge *edg = elm_side->edge();

			for(int j=0; j < edg->n_sides;j++) {
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					unsigned int sousedni_element = other_side->element()->index(); // 2D element

							//xprintf(Msg, "Naleznut sousedni element elementu(3030) s ctyrstenem(%d)\n", ele->index());
							//xprintf(Msg, "\t\t Idx původního elementu a jeho hrany(%d,%d) - Idx sousedního elementu a jeho hrany(%d,%d)\n",elm->index(),stena,other_side->element()->index(),other_side->el_idx());


						if(!intersectionExists(sousedni_element,elm_2D)){
							//flag_for_3D_elements[ele->index()] = sousedni_element;

							// Vytvoření průniku bez potřeby počítání
							IntersectionLocal il_other(sousedni_element, elm_2D);
							intersection_list[sousedni_element].push_back(il_other);

							ProlongationLine pl2(sousedni_element, elm_2D, intersection_list[sousedni_element].size() - 1);
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
					//

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

	//xprintf(Msg, "elementy(%d,%d)\n", pl.getElement2DIdx(), pl.getElement3DIdx());

	update_triangle(elm);
	update_tetrahedron(ele);

	/*if(elm->index() == 5648 && ele->index() == 180849){
		triangle.to_string();
		tetrahedron.to_string();
		//return;
	}*/

	/*if(elm->index() == 1076 && ele->index() == 223421){

		triangle.to_string();
		tetrahedron.to_string();

	}

	xprintf(Msg, "TU(%d,%d)\n",elm->index(),ele->index());*/

	element_2D_index = pl.getElement2DIdx();

	/*if(element_2D_index == 3030){
		xprintf(Msg, "OU YEAH - 3030 a 3D(%d)\n",pl.getElement3DIdx());
		if(pl.getElement3DIdx() == 36762){
			xprintf(Msg, "no sakra\n");
		}
	}*/

	ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
	CI_23.init();
	CI_23.compute(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()]);

	/*if(elm->index() == 45 && ele->index() == 5874){
		xprintf(Msg,"ou yeah\n");
		//intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].print();
	}*/


	if(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].size() > 2){


		/*if(elm->index() == 45 && ele->index() == 5874){
				xprintf(Msg,"ou yeah 2\n");
				xprintf(Msg,"Velikost %d\n",intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].size());
				xprintf(Msg,"Patologicky: %d\n",intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].isPatological());
				intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].print();
		}*/


		std::vector<unsigned int> prolongation_table;
		intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()].traceGenericPolygon(prolongation_table);
		/*if(elm->index() == 45 && ele->index() == 5874){
			xprintf(Msg,"neco \n");
			for(unsigned int aa = 0; aa < prolongation_table.size();aa++){
				xprintf(Msg, "PT %d : %d\n", aa, prolongation_table[aa]);
			}
		}*/
		//all_intersections.push_back(intersection_list[pl.getElement2DIdx()][pl.getDictionaryIdx()]);
		computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);
		/*if(elm->index() == 1076 && ele->index() == 223421){
						//xprintf(Msg,"ou yeah 4\n");
							}*/
	}
};


void InspectElements::ComputeIntersections23(){

	{ START_TIMER("Incializace pruniku");
	compute_intersections_init<2,3>();
	END_TIMER("Inicializace pruniku");}

	{ START_TIMER("Prochazeni vsech elementu");
	FOR_ELEMENTS(mesh, elm) {
		if (elm->dim() == 2 && !closed_elements[elm.index()] && elements_bb[elm->index()].intersect(mesh_3D_bb)) {
			update_triangle(elm);

			{ START_TIMER("Hlavni vypocet");

			FOR_ELEMENTS(mesh, ele) {
				if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] != (int)(elm->index()) && elements_bb[elm->index()].intersect(elements_bb[ele->index()])) {
					update_tetrahedron(ele);
					IntersectionLocal il(elm.index(), ele->index());
					ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
					CI_23.init();
					CI_23.compute(il);

					if(il.size() > 2){

						std::vector<unsigned int> prolongation_table;
						il.traceGenericPolygon(prolongation_table);

						{ START_TIMER("Prochazeni vsech front");

						intersection_list[elm.index()].push_back(il);
						flag_for_3D_elements[ele->index()] = elm.index();
						//xprintf(Msg,"Velikost intersection_list(%d), pocet IL v konkretnim listu(%d)\n", intersection_list.size(), intersection_list[elm.index()].size());

						// PRODLUŽOVÁNÍ
						computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);

						while(1){
							while(!prolongation_line_queue_3D.empty()){
								computeIntersections2d3dProlongation((ProlongationLine)prolongation_line_queue_3D.front());
								prolongation_line_queue_3D.pop();
							}

							if(element_2D_index >= 0){
								closed_elements[element_2D_index] = true;
							}

							element_2D_index = -1;
							// Pridat priznak trojuhleniku, ze je projity
							if(!prolongation_line_queue_2D.empty()){

								computeIntersections2d3dProlongation((ProlongationLine)prolongation_line_queue_2D.front());
								prolongation_line_queue_2D.pop();

							}

							if(element_2D_index >= 0){
								closed_elements[element_2D_index] = true;

							}

							element_2D_index = -1;

							if(prolongation_line_queue_2D.empty() && prolongation_line_queue_3D.empty()){

								break;
							}

						}
						END_TIMER("Prochazeni vsech front");}
						//prunik = true;
						break; // ukončí procházení dalších bounding boxů
					}

				}
			}
			// Prošlo se celé pole sousedním bounding boxů, pokud nevznikl průnik, může se trojúhelník nastavit jako uzavřený
			closed_elements[elm.index()] = true;
			END_TIMER("Hlavni vypocet");}

		}
	}

	END_TIMER("Prochazeni vsech elementu");}

};


bool InspectElements::intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx){

	bool found = false;

	for(unsigned int i = 0; i < intersection_list[elm_2D_idx].size();i++){

		if(intersection_list[elm_2D_idx][i].idx_3D() == elm_3D_idx){
			found = true;
			break;
		}

	}
	return found;
};

void InspectElements::update_abscissa(const ElementFullIter &element_1D){
	arma::vec3 *field_of_points[2] = {&element_1D->node[0]->point(),&element_1D->node[1]->point()};
	abscissa.setSimplices(field_of_points);
};

void InspectElements::update_triangle(const ElementFullIter &element_2D){
	arma::vec3 *field_of_points[3] = {&element_2D->node[0]->point(),&element_2D->node[1]->point(),&element_2D->node[2]->point()};
	triangle.setSimplices(field_of_points);
};

void InspectElements::update_tetrahedron(const ElementFullIter &element_3D){
	arma::vec3 *field_of_points[4] = {&element_3D->node[0]->point(),&element_3D->node[1]->point(),&element_3D->node[2]->point(),&element_3D->node[3]->point()};
	tetrahedron.setSimplices(field_of_points);
};

void InspectElements::print_mesh_to_file(string name){


	for(unsigned int i = 0; i < 2;i++){
		string t_name = name;

		unsigned int number_of_intersection_points = 0;
		unsigned int number_of_polygons = 0;
		unsigned int number_of_nodes = mesh->n_nodes();

		for(unsigned int j = 0; j < intersection_list.size();j++){
			number_of_polygons += intersection_list[j].size();
			for(unsigned int k = 0; k < intersection_list[j].size();k++){
				number_of_intersection_points += intersection_list[j][k].size();
			}
		}

		FILE * file;
		if(i == 0){
			file = fopen((t_name.append("_0.msh")).c_str(),"w");
		}else{
			file = fopen((t_name.append("_1.msh")).c_str(),"w");
		}
		fprintf(file, "$MeshFormat\n");
		fprintf(file, "2.2 0 8\n");
		fprintf(file, "$EndMeshFormat\n");
		fprintf(file, "$Nodes\n");
		fprintf(file, "%d\n", (mesh->n_nodes() + number_of_intersection_points));

		FOR_NODES(mesh, nod){
			arma::vec3 _nod = nod->point();
			fprintf(file,"%d %f %f %f\n", nod.id(), _nod[0], _nod[1], _nod[2]);
		}

		for(unsigned int j = 0; j < intersection_list.size();j++){
			for(unsigned int k = 0; k < intersection_list[j].size();k++){

				IntersectionLocal il = intersection_list[j][k];

				ElementFullIter el2D = mesh->element(il.idx_2D());
				ElementFullIter el3D = mesh->element(il.idx_3D());


				for(unsigned int l = 0; l < il.size(); l++){
					//xprintf(Msg, "first\n");
					number_of_nodes++;
					IntersectionPoint<2,3> IP23 = il.get_point(l);
					arma::vec3 _global;
					if(i == 0){
						_global = (IP23.get_local_coords1())[0] * el2D->node[0]->point()
														   +(IP23.get_local_coords1())[1] * el2D->node[1]->point()
														   +(IP23.get_local_coords1())[2] * el2D->node[2]->point();
					}else{
						_global = (IP23.get_local_coords2())[0] * el3D->node[0]->point()
														   +(IP23.get_local_coords2())[1] * el3D->node[1]->point()
														   +(IP23.get_local_coords2())[2] * el3D->node[2]->point()
														   +(IP23.get_local_coords2())[3] * el3D->node[3]->point();
					}
					fprintf(file,"%d %f %f %f\n", number_of_nodes, _global[0], _global[1], _global[2]);
				}
			}
		}

		fprintf(file,"$EndNodes\n");
		fprintf(file,"$Elements\n");
		fprintf(file,"%d\n", (number_of_intersection_points + mesh->n_elements()) );

		FOR_ELEMENTS(mesh, elee){
			// 1 4 2 30 26 1 2 3 4
			// 2 2 2 2 36 5 6 7
			if(elee->dim() == 3){
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				int id3 = mesh->node_vector.index(elee->node[2]) + 1;
				int id4 = mesh->node_vector.index(elee->node[3]) + 1;

				fprintf(file,"%d 4 2 30 26 %d %d %d %d\n", elee.id(), id1, id2, id3, id4);
			}else if(elee->dim() == 2){
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				int id3 = mesh->node_vector.index(elee->node[2]) + 1;
				fprintf(file,"%d 2 2 2 36 %d %d %d\n", elee.id(), id1, id2, id3);

			}else{
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				fprintf(file,"%d 1 2 14 16 %d %d\n",elee.id(), id1, id2);
				//xprintf(Msg, "Missing implementation of printing 1D element to a file\n");
			}
		}


		unsigned int number_of_elements = mesh->n_elements();
		unsigned int last = 0;
		unsigned int nodes = mesh->n_nodes();

		for(unsigned int j = 0; j < intersection_list.size();j++){
				for(unsigned int k = 0; k < intersection_list[j].size();k++){

				IntersectionLocal il = intersection_list[j][k];

				for(unsigned int l = 0; l < il.size(); l++){
					number_of_elements++;
					nodes++;
					if(l == 0){
						last = nodes;
					}

					if((l+1) == il.size()){
						fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, last);
					}else{
						fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes+1);
					}

				}
			}
		}

		fprintf(file,"$EndElements\n");
		fclose(file);
	}

};

void InspectElements::print_mesh_to_file_1D(string name){


	for(unsigned int i = 0; i < 2;i++){
		string t_name = name;

		unsigned int number_of_intersection_points = 0;
		unsigned int number_of_polygons = 0;
		unsigned int number_of_nodes = mesh->n_nodes();
		//xprintf(Msg, "Zde6\n");
		for(unsigned int j = 0; j < intersection_line_list.size();j++){
			number_of_polygons += intersection_line_list[j].size();
			for(unsigned int k = 0; k < intersection_line_list[j].size();k++){
				number_of_intersection_points += intersection_line_list[j][k].size();
			}
		}
		//xprintf(Msg, "Zde5\n");
		FILE * file;
		if(i == 0){
			file = fopen((t_name.append("_0.msh")).c_str(),"w");
		}else{
			file = fopen((t_name.append("_1.msh")).c_str(),"w");
		}
		fprintf(file, "$MeshFormat\n");
		fprintf(file, "2.2 0 8\n");
		fprintf(file, "$EndMeshFormat\n");
		fprintf(file, "$Nodes\n");
		fprintf(file, "%d\n", (mesh->n_nodes() + number_of_intersection_points));
		//xprintf(Msg, "Zde4\n");
		FOR_NODES(mesh, nod){
			arma::vec3 _nod = nod->point();
			fprintf(file,"%d %f %f %f\n", nod.id(), _nod[0], _nod[1], _nod[2]);
		}
		//xprintf(Msg, "Zde3\n");
		for(unsigned int j = 0; j < intersection_line_list.size();j++){
			for(unsigned int k = 0; k < intersection_line_list[j].size();k++){
				//xprintf(Msg, "Zde3.1\n");
				IntersectionLine il = intersection_line_list[j][k];
				//xprintf(Msg, "Zde3.2\n");
				//xprintf(Msg, "el %d %d\n", il.get_elm1D_idx(), il.get_elm3D_idx());
				ElementFullIter el1D = mesh->element(il.get_elm1D_idx());
				ElementFullIter el3D = mesh->element(il.get_elm3D_idx());
				//xprintf(Msg, "Zde3.3\n");

				for(unsigned int l = 0; l < il.size(); l++){
					//xprintf(Msg, "Zde3.4\n");
					//xprintf(Msg, "first\n");
					number_of_nodes++;
					IntersectionPoint<1,3> IP13 = il.get_point(l);
					arma::vec3 _global;
					if(i == 0){
						_global = (IP13.get_local_coords1())[0] * el1D->node[0]->point()
														   +(IP13.get_local_coords1())[1] * el1D->node[1]->point();
					}else{
						_global = (IP13.get_local_coords2())[0] * el3D->node[0]->point()
														   +(IP13.get_local_coords2())[1] * el3D->node[1]->point()
														   +(IP13.get_local_coords2())[2] * el3D->node[2]->point()
														   +(IP13.get_local_coords2())[3] * el3D->node[3]->point();
					}
					//xprintf(Msg, "Zde3.5\n");
					fprintf(file,"%d %f %f %f\n", number_of_nodes, _global[0], _global[1], _global[2]);
				}
			}
		}

		fprintf(file,"$EndNodes\n");
		fprintf(file,"$Elements\n");
		fprintf(file,"%d\n", (number_of_intersection_points/2 + mesh->n_elements()) );
		//xprintf(Msg, "Zde2\n");
		FOR_ELEMENTS(mesh, elee){
			// 1 4 2 30 26 1 2 3 4
			// 2 2 2 2 36 5 6 7
			if(elee->dim() == 3){
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				int id3 = mesh->node_vector.index(elee->node[2]) + 1;
				int id4 = mesh->node_vector.index(elee->node[3]) + 1;

				fprintf(file,"%d 4 2 30 26 %d %d %d %d\n", elee.id(), id1, id2, id3, id4);
			}else if(elee->dim() == 2){
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				int id3 = mesh->node_vector.index(elee->node[2]) + 1;
				fprintf(file,"%d 2 2 2 36 %d %d %d\n", elee.id(), id1, id2, id3);

			}else if(elee->dim() == 1){
				int id1 = mesh->node_vector.index(elee->node[0]) + 1;
				int id2 = mesh->node_vector.index(elee->node[1]) + 1;
				fprintf(file,"%d 1 2 14 16 %d %d\n",elee.id(), id1, id2);
				//xprintf(Msg, "Missing implementation of printing 1D element to a file\n");
			}
		}


		//xprintf(Msg, "Zde\n");

		unsigned int number_of_elements = mesh->n_elements();
		unsigned int last = 0;
		unsigned int nodes = mesh->n_nodes();

		for(unsigned int j = 0; j < intersection_line_list.size();j++){
				for(unsigned int k = 0; k < intersection_line_list[j].size();k++){

				IntersectionLine il = intersection_line_list[j][k];


				if(il.size() == 1){
					number_of_elements++;
					nodes++;
					fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes);

				}else if(il.size() == 2){
					number_of_elements++;
					nodes++;
					fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes+1);
					nodes++;
				}


			}
		}

		fprintf(file,"$EndElements\n");
		fclose(file);
	}

};

} // END namespace
