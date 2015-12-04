/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "inspectelements.h"
#include "prolongation.h"
#include "intersectionpoint.h"
#include "intersectionline.h"
#include "intersectionpolygon.h"
#include "computeintersection.h"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/bih_tree.hh"

#include "mesh/ngh/include/triangle.h"
#include "mesh/ngh/include/abscissa.h"

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
		elements_bb.resize(mesh->n_elements());
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
	intersection_list.assign(mesh->n_elements(),std::vector<IntersectionPolygon>());

	if(elements_bb.size() == 0){
		elements_bb.resize(mesh->n_elements());
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
        if (elm->dim() == 1 &&                                  // is 1D element
            !closed_elements[elm->index()] &&                   // is not closed yet
            elements_bb[elm->index()].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {
            DBGMSG("-----Nalezen 1D element------ \n");

			update_abscissa(elm);
			std::vector<unsigned int> searchedElements;
			bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

            // go through all 3D elements that can have possibly intersection with 1D elements bounding box
			for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
				int idx = *it;
				ElementFullIter ele = mesh->element( idx );
				if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] == -1) {

					update_tetrahedron(ele);

					IntersectionLine il(elm->index(), ele->index());
					ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
					CI_13.init();
					CI_13.compute(il.points());

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
	for(const IntersectionPoint<1,3> &IP : il.points()) {
        std::cout << IP;        
		if(IP.dim_A() == 0) { 
            // if IP is end of the 1D element
            // prolongate 1D element as long as it creates prolongation point on the side of tetrahedron

            SideIter elm_side = elm->side(IP.idx_A());  // side of 1D element is vertex
			Edge *edg = elm_side->edge();
			for(int j=0; j < edg->n_sides;j++) {

				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {

					unsigned int sousedni_element = other_side->element()->index();
					if(!closed_elements[sousedni_element]){
						prolongate_1D_element(other_side->element(), ele);  //computes intersection with the same tetrahedron
					}
				}
			}
		}else{
            std::vector<Edge*> edges;

            switch (IP.dim_B())
            {
                // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
                case 0: for(unsigned int j=0; j < RefElement<3>::n_sides_per_node; j++)
                            edges.push_back(&(mesh->edges[ele->edge_idx_[RefElement<3>::interact<2,0>(IP.idx_B())[j]]]));
                        break;
                
                // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
                case 1: for(unsigned int j=0; j < RefElement<3>::n_sides_per_line; j++)
                            edges.push_back(&(mesh->edges[ele->edge_idx_[RefElement<3>::interact<2,1>(IP.idx_B())[j]]]));
                        break;
                        
                // IP is on a side of tetrahedron; only possible edge is from the given side
                case 2: edges.push_back(&(mesh->edges[ele->edge_idx_[IP.idx_B()]]));
                        break;
                default: ASSERT_LESS(IP.dim_B(),3);
            }
            
            for(Edge* edg : edges)
            for(int j=0; j < edg->n_sides;j++) {
                if (edg->side(j)->element() != ele) {
                    unsigned int sousedni_element = edg->side(j)->element()->index();

                    if(flag_for_3D_elements[sousedni_element] == -1){
                        flag_for_3D_elements[sousedni_element] = elm->index();
                        ProlongationPoint pp = {elm->index(), sousedni_element, ele->index()};
                        prolongation_point_queue.push(pp);
                    }
                }   
            }
        }  
//         if(il[i].is_vertex()){
// 
//             SideIter elm_side = elm->side((unsigned int)(1-il[i].local_bcoords_A()[0]));
//             Edge *edg = elm_side->edge();
//             for(int j=0; j < edg->n_sides;j++) {
// 
//                 SideIter other_side=edg->side(j);
//                 if (other_side != elm_side) {
// 
//                     unsigned int sousedni_element = other_side->element()->index();
//                     if(!closed_elements[sousedni_element]){
//                         prolongate_1D_element(other_side->element(), ele);
//                     }
//                 }
//             }
//         }else{
// 			SideIter elm_side = ele->side((unsigned int)(3-il[i].idx_B())); //ele->side(3-stena);
// 			
// 			Edge *edg = elm_side->edge();
//             
// 			for(int j=0; j < edg->n_sides;j++) {
// 
// 				SideIter other_side=edg->side(j);
// 				if (other_side != elm_side) {
// 					unsigned int sousedni_element = other_side->element()->index();
// 
// 					if(flag_for_3D_elements[sousedni_element] == -1){ // || (flag_for_3D_elements[sousedni_element] != (int)elm->index() && !intersectionExists(elm->index(),sousedni_element))){
// 						flag_for_3D_elements[sousedni_element] = elm->index();
// 						ProlongationPoint pp = {elm->index(), sousedni_element, ele->index()};
// 						prolongation_point_queue.push(pp);
// 					}
// 				}
// 			}
// 		}
	}
};

void InspectElements::prolongate_1D_element(const ElementFullIter &elm, const ElementFullIter &ele){

    update_abscissa(elm);
    closed_elements[elm->index()] = true;

    IntersectionLine il(elm->index(), ele->index());
    ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
    CI_13.init();
    CI_13.compute(il.points());

    if(il.size() > 1){
        intersection_line_list[elm->index()].push_back(il);
        prolongate_elements(il, elm, ele);
    }
};

void InspectElements::prolongate(const ProlongationPoint &pp){

    ElementFullIter elm = mesh->element(pp.elm_1D_idx);
    ElementFullIter ele = mesh->element(pp.elm_3D_idx);

    update_abscissa(elm);
    update_tetrahedron(ele);

    IntersectionLine il(elm->index(), ele->index());
    ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
    CI_13.init();
    CI_13.compute(il.points());

    if(il.size() > 1){
        intersection_line_list[elm->index()].push_back(il);
        prolongate_elements(il, elm, ele);
    }

}

template<>
void InspectElements::compute_intersections<2,3>(){
	{ START_TIMER("Incializace pruniku");
		compute_intersections_init<2,3>();
		END_TIMER("Inicializace pruniku");}

		BIHTree bt(mesh, 20);
		//bool nalezen = false;

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
						IntersectionPolygon il(elm.index(), ele->index());
						ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
						CI_23.init();
						CI_23.compute(il);

						if(il.size() > 2){

							//nalezen = true;
							std::vector<unsigned int> prolongation_table;
							il.trace_polygon(prolongation_table);
							//il.printTracingTable();

							{ START_TIMER("Prochazeni vsech front");

							intersection_list[elm.index()].push_back(il);
							flag_for_3D_elements[ele->index()] = elm.index();
                            DBGMSG("Velikost intersection_list(%d), pocet IL v konkretnim listu(%d)\n", 
                                   intersection_list.size(), intersection_list[elm.index()].size());

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
				/*if(nalezen){
					break;
				}*/
				END_TIMER("Hlavni vypocet");}

			}
		}

		END_TIMER("Prochazeni vsech elementu");}

};


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
							IntersectionPolygon il_other(sousedni_element, elm_2D);
							intersection_list[sousedni_element].push_back(il_other);

							//ProlongationLine pl2(sousedni_element, elm_2D, intersection_list[sousedni_element].size() - 1);
							ProlongationLine pl2 = {sousedni_element, elm_2D, intersection_list[sousedni_element].size() - 1, -1, -1};
							prolongation_line_queue_2D.push(pl2);
						}
				}
			}

		}else{
			// prodlužuji stěnou

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
						IntersectionPolygon il_other(elm.index(), sousedni_element);
						intersection_list[elm.index()].push_back(il_other);

						//ProlongationLine pl2(elm.index(), sousedni_element, intersection_list[elm.index()].size() - 1);
						ProlongationLine pl2 = {elm.index(), sousedni_element, intersection_list[elm.index()].size() - 1, -1, -1};
						prolongation_line_queue_3D.push(pl2);

					}
				}
			}

		}

	}
};

void InspectElements::computeIntersections2d3dProlongation(const ProlongationLine &pl){

	// Měly bychom být vždy ve stejném trojúhelníku, tudíž updateTriangle by měl být zbytečný
	ElementFullIter elm = mesh->element(pl.elm_2D_idx);
	ElementFullIter ele = mesh->element(pl.elm_3D_idx);

	update_triangle(elm);
	update_tetrahedron(ele);


	element_2D_index = pl.elm_2D_idx;


	ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
	CI_23.init();
	CI_23.compute(intersection_list[pl.elm_2D_idx][pl.dictionary_idx]);



	if(intersection_list[pl.elm_2D_idx][pl.dictionary_idx].size() > 2){

		std::vector<unsigned int> prolongation_table;
		intersection_list[pl.elm_2D_idx][pl.dictionary_idx].trace_polygon(prolongation_table);
		computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);
	}
};


bool InspectElements::intersectionExists(unsigned int elm_2D_idx, unsigned int elm_3D_idx){

	bool found = false;

	for(unsigned int i = 0; i < intersection_list[elm_2D_idx].size();i++){

		if(intersection_list[elm_2D_idx][i].ele_3d_idx() == elm_3D_idx){
			found = true;
			break;
		}

	}
	return found;
};

void InspectElements::update_abscissa(const ElementFullIter &element_1D){
	arma::vec3 *field_of_points[2] = {&element_1D->node[0]->point(),&element_1D->node[1]->point()};
	abscissa.set_simplices(field_of_points);
};

void InspectElements::update_triangle(const ElementFullIter &element_2D){
	arma::vec3 *field_of_points[3] = {&element_2D->node[0]->point(),&element_2D->node[1]->point(),&element_2D->node[2]->point()};
	triangle.set_simplices(field_of_points);
};

void InspectElements::update_tetrahedron(const ElementFullIter &element_3D){
	arma::vec3 *field_of_points[4] = {&element_3D->node[0]->point(),&element_3D->node[1]->point(),&element_3D->node[2]->point(),&element_3D->node[3]->point()};
	tetrahedron.set_simplices(field_of_points);
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

				IntersectionPolygon il = intersection_list[j][k];

				ElementFullIter el2D = mesh->element(il.ele_2d_idx());
				ElementFullIter el3D = mesh->element(il.ele_3d_idx());


				for(unsigned int l = 0; l < il.size(); l++){
					//xprintf(Msg, "first\n");
					number_of_nodes++;
					IntersectionPoint<2,3> IP23 = il[l];
					arma::vec3 _global;
					if(i == 0){
						_global = (IP23.local_bcoords_A())[0] * el2D->node[0]->point()
														   +(IP23.local_bcoords_A())[1] * el2D->node[1]->point()
														   +(IP23.local_bcoords_A())[2] * el2D->node[2]->point();
					}else{
						_global = (IP23.local_bcoords_B())[0] * el3D->node[0]->point()
														   +(IP23.local_bcoords_B())[1] * el3D->node[1]->point()
														   +(IP23.local_bcoords_B())[2] * el3D->node[2]->point()
														   +(IP23.local_bcoords_B())[3] * el3D->node[3]->point();
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

				IntersectionPolygon il = intersection_list[j][k];

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
				ElementFullIter el1D = mesh->element(il.ele_1d_idx());
				ElementFullIter el3D = mesh->element(il.ele_3d_idx());
				//xprintf(Msg, "Zde3.3\n");

				for(unsigned int l = 0; l < il.size(); l++){
					//xprintf(Msg, "Zde3.4\n");
					//xprintf(Msg, "first\n");
					number_of_nodes++;
					IntersectionPoint<1,3> IP13 = il[l];
					arma::vec3 _global;
					if(i == 0){
						_global = (IP13.local_bcoords_A())[0] * el1D->node[0]->point()
														   +(IP13.local_bcoords_A())[1] * el1D->node[1]->point();
					}else{
						_global = (IP13.local_bcoords_B())[0] * el3D->node[0]->point()
														   +(IP13.local_bcoords_B())[1] * el3D->node[1]->point()
														   +(IP13.local_bcoords_B())[2] * el3D->node[2]->point()
														   +(IP13.local_bcoords_B())[3] * el3D->node[3]->point();
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
		//unsigned int last = 0;
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

double InspectElements::polygonArea()
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_list.size(); i++){
        for(unsigned int j = 0; j < intersection_list[i].size();j++){
            Element efi = *mesh->element(intersection_list[i][j].ele_2d_idx());
             TTriangle t2d(efi);
             double t2dArea = t2d.GetArea();
             double localArea = intersection_list[i][j].get_area();//il.getArea();
             subtotal += 2*localArea*t2dArea;
        }
    }
    return subtotal;
}

double InspectElements::line_length()
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_line_list.size(); i++){
        for(unsigned int j = 0; j < intersection_line_list[i].size();j++){
            Element efi = *mesh->element(intersection_line_list[i][j].ele_1d_idx());
            TAbscissa t1d(efi);
            double t1d_length = t1d.Length();
            DBGMSG("t1d length %f\n",t1d_length);
            double local_length = intersection_line_list[i][j].compute_length();
            subtotal += local_length*t1d_length;
        }
    }
    return subtotal;
}


} // END namespace
