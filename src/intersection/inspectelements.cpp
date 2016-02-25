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

namespace computeintersection {

InspectElements::InspectElements(Mesh* _mesh):mesh(_mesh){
	element_2D_index = -1;
};

InspectElements::~InspectElements(){};

template<>
void InspectElements::compute_intersections_init<1,2>(){
    flag_for_3D_elements.assign(mesh->n_elements(), -1);
    closed_elements.assign(mesh->n_elements(), false);
    intersection_point_list.assign(mesh->n_elements(),std::vector<IntersectionPoint<1,2>>());
    
    if(elements_bb.size() == 0){
        elements_bb.resize(mesh->n_elements());
        bool first_2d_element = true;
        FOR_ELEMENTS(mesh, elm) {

            elements_bb[elm->index()] = elm->bounding_box();

                if (elm->dim() == 2){
                    if(first_2d_element){
                        first_2d_element = false;
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
void InspectElements::compute_intersections<1,2>(){

    compute_intersections_init<1,2>();
    BIHTree bt(mesh, 20);

    FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == 1 &&                                  // is 1D element
            !closed_elements[elm->index()] &&                   // is not closed yet
            elements_bb[elm->index()].intersect(mesh_3D_bb))    // its bounding box intersects 2D mesh bounding box
        {
            DBGMSG("-----Nalezen 1D element------ \n");

            update_abscissa(elm);
            std::vector<unsigned int> searchedElements;
            bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

            // go through all 2D elements that can have possibly intersection with 1D elements bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
                int idx = *it;
                ElementFullIter ele = mesh->element( idx );
                if (ele->dim() == 2 && flag_for_3D_elements[ele->index()] == -1) {

                    update_triangle(ele);

//                     vector<IntersectionPoint> IPs;
//                     IPs.reserve(2);
                    START_TIMER("Compute intersection 12");
                    ComputeIntersection<Simplex<1>, Simplex<2>> CI_12(abscissa, triangle);
                    bool ip_found = CI_12.compute_final(intersection_point_list[elm->index()]);
                    END_TIMER("Compute intersection 12");
                    
//                     DBGMSG("IL IPs:\n");
//                     for(unsigned int i=0; i < il.size(); i++)
//                         std::cout << il.points()[i];
                    
                    if(ip_found){
                        closed_elements[elm->index()] = true;
                        flag_for_3D_elements[ele->index()] = elm->index();

                        //TODO: prolongation in parallel pathologic case (2 IPs) (I dont know what about IP being triangle vertex)
//                         prolongate_elements(il, elm, ele);
// 
//                         while(!prolongation_point_queue.empty()){
// 
//                             prolongate((ProlongationPoint)prolongation_point_queue.front());
//                             prolongation_point_queue.pop();
// 
//                         }
                        break;
                    }
                }
            }
        }
    }
}

template<>
void InspectElements::compute_intersections<1,3>(){

	compute_intersections_init<1,3>();
	BIHTree bt(mesh, 20);

    // Find the first candidates (1D elements) for intersection.
	FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == 1 &&                                  // is 1D element
            !closed_elements[elm->index()] &&                   // is not closed yet
            elements_bb[elm->index()].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {
            DBGMSG("----- Found 1D element candidate for intersection [idx = %d]------ \n",elm->index());
            flag_for_3D_elements.assign(mesh->n_elements(), -1);
            
			update_abscissa(elm);
            // create vector of elements which bounding box interact with bounding box of 1D element
			std::vector<unsigned int> searchedElements;
			bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

            // go through all 3D elements that can have possibly intersection with 1D elements bounding box
			for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++){
				int idx = *it;
				ElementFullIter ele = mesh->element( idx );
				if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] == -1) {

					update_tetrahedron(ele);

					IntersectionLine il(elm->index(), ele->index());
                    START_TIMER("Compute intersection 13");
					ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
					CI_13.init();
					CI_13.compute(il.points());
                    END_TIMER("Compute intersection 13");

//                     DBGMSG("IL IPs:\n");
//                     for(unsigned int i=0; i < il.size(); i++)
//                         std::cout << il.points()[i];
                    
					if(il.size() > 1){
						closed_elements[elm->index()] = true;
						flag_for_3D_elements[ele->index()] = elm->index();
						intersection_line_list[elm->index()].push_back(il);

						prolongate_1D_decide(il, elm, ele);

						while(!prolongation_point_queue.empty()){
                            ProlongationPoint pp = prolongation_point_queue.front();
                            DBGMSG("Prolongation queue compute in 3d ele_idx %d.\n",pp.elm_3D_idx);
                            
                            ElementFullIter e1d = mesh->element(pp.elm_1D_idx),
                                            e3d = mesh->element(pp.elm_3D_idx);
                            update_abscissa(e1d);
                            update_tetrahedron(e3d);
                            prolongate_1D_element(e1d, e3d);
                            
							prolongation_point_queue.pop();
						}
						break;
					}
				}
			}
		}
	}
}

void InspectElements::prolongate_1D_decide(const computeintersection::IntersectionLine& il, const ElementFullIter& ele_1d, const ElementFullIter& ele_3d){
	for(const IntersectionPoint<1,3> &IP : il.points()) {
		if(IP.dim_A() == 0) { 
            // if IP is end of the 1D element
            // prolongate 1D element as long as it creates prolongation point on the side of tetrahedron
            SideIter elm_side = ele_1d->side(IP.idx_A());  // side of 1D element is vertex
			Edge *edg = elm_side->edge();
			for(int j=0; j < edg->n_sides;j++) {
                
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {

					unsigned int sousedni_element = other_side->element()->index();
					if(!closed_elements[sousedni_element]){
                        DBGMSG("Call prolongate 1D to %d.\n", sousedni_element);
                        update_abscissa(other_side->element());
						prolongate_1D_element(other_side->element(), ele_3d);  //computes intersection with the same tetrahedron
					}
				}
			}
		}else{
            std::vector<Edge*> edges;

            switch (IP.dim_B())
            {
                // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
                case 0: for(unsigned int j=0; j < RefElement<3>::n_sides_per_node; j++)
                            edges.push_back(&(mesh->edges[ele_3d->edge_idx_[RefElement<3>::interact<2,0>(IP.idx_B())[j]]]));
                        break;
                
                // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
                case 1: for(unsigned int j=0; j < RefElement<3>::n_sides_per_line; j++)
                            edges.push_back(&(mesh->edges[ele_3d->edge_idx_[RefElement<3>::interact<2,1>(IP.idx_B())[j]]]));
                        break;
                        
                // IP is on a side of tetrahedron; only possible edge is from the given side
                case 2: edges.push_back(&(mesh->edges[ele_3d->edge_idx_[IP.idx_B()]]));
                        break;
                default: ASSERT_LESS(IP.dim_B(),3);
            }
            
            for(Edge* edg : edges)
            for(int j=0; j < edg->n_sides;j++) {
                if (edg->side(j)->element() != ele_3d) {
                    unsigned int sousedni_element = edg->side(j)->element()->index();

                    if(flag_for_3D_elements[sousedni_element] == -1){
                        DBGMSG("Prolongate push to queue %d.\n", sousedni_element);
                        flag_for_3D_elements[sousedni_element] = ele_1d->index();
                        ProlongationPoint pp = {ele_1d->index(), sousedni_element, ele_3d->index()};
                        prolongation_point_queue.push(pp);
                    }
                }   
            }
        }  
	}
};

void InspectElements::prolongate_1D_element(const ElementFullIter& ele_1d, const ElementFullIter& ele_3d){

    DBGMSG("Prolongate 1D: %d in %d.\n", ele_1d->index(), ele_3d->index());
    closed_elements[ele_1d->index()] = true;

    IntersectionLine il(ele_1d->index(), ele_3d->index());
    ComputeIntersection<Simplex<1>, Simplex<3>> CI_13(abscissa, tetrahedron);
    CI_13.init();
    CI_13.compute(il.points());

    if(il.size() > 1){
        //cout << il[0] << il[1];
        intersection_line_list[ele_1d->index()].push_back(il);
        prolongate_1D_decide(il, ele_1d, ele_3d);
    }
};

template<>
void InspectElements::compute_intersections<2,3>(){
	{   START_TIMER("Intersection initialization");
		compute_intersections_init<2,3>();
		END_TIMER("Intersection initialization");}

		BIHTree bt(mesh, 20);
		//bool nalezen = false;

		{ START_TIMER("Element iteration");
		FOR_ELEMENTS(mesh, elm) {
			if (elm->dim() == 2 &&                                  // is 2D element
                !closed_elements[elm.index()] &&                    // is not closed yet
                elements_bb[elm->index()].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
            {    
				update_triangle(elm);
				std::vector<unsigned int> searchedElements;
				bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

				{ START_TIMER("Bounding box element iteration");
				for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
				{
					unsigned int idx = *it;
					ElementFullIter ele = mesh->element(idx);

					if (ele->dim() == 3 && flag_for_3D_elements[ele->index()] != (int)(elm->index())) {
						update_tetrahedron(ele);
						IntersectionPolygon il(elm.index(), ele->index());
                        START_TIMER("Compute intersection 23");
						ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
						CI_23.init();
						CI_23.compute(il);
                        END_TIMER("Compute intersection 23");

                        //TODO: how does prolongation work when there are 1 or 2 IPs?
						if(il.size() > 2){

							//nalezen = true;
							std::vector<unsigned int> prolongation_table;
							il.trace_polygon(prolongation_table);
							//il.printTracingTable();

							{ START_TIMER("Prolongation queue");

							intersection_list[elm.index()].push_back(il);
							flag_for_3D_elements[ele->index()] = elm.index();
                            DBGMSG("Velikost intersection_list(%d), pocet IL v konkretnim listu(%d)\n", 
                                   intersection_list.size(), intersection_list[elm.index()].size());

							// PRODLUŽOVÁNÍ
							computeIntersections2d3dUseProlongationTable(prolongation_table, elm, ele);

                            // - find first intersection polygon
                            // - fill both prolongation queue in 2D and 3D
                            // - clear prolongation queue:
                            //      - clear 3D prolongation queue:
                            //          - compute CI
                            //          - fill both prolongation queue 2D and 3D
                            //          - repeat until 3D queue is empty
                            //      - we can close the 2D element which we started with
                            //      - clear 2D prolongation queue:
                            //          - compute CI
                            //          - fill both prolongation queue 2D and 3D
                            //          - repeat until 2D queue is empty
                            //      - emptying 2D queue could fill 3D queue again, so repeat
							while(1){
								while(!prolongation_line_queue_3D.empty()){
                                    ProlongationLine pl = prolongation_line_queue_3D.front();
                                    DBGMSG("Prolongation queue compute in 3d ele_idx %d.\n",pl.elm_3D_idx);
									computeIntersections2d3dProlongation(pl);
									prolongation_line_queue_3D.pop();

								}

								if(element_2D_index >= 0){
									closed_elements[element_2D_index] = true;
								}

								element_2D_index = -1;
								if(!prolongation_line_queue_2D.empty()){
                                    ProlongationLine pl = prolongation_line_queue_2D.front();
                                    DBGMSG("Prolongation queue compute in 2d ele_idx %d.\n",pl.elm_2D_idx);
									computeIntersections2d3dProlongation(pl);
									prolongation_line_queue_2D.pop();

								}

								// TODO shoud not this be further, just before 'break' when queues are empty??
								if(element_2D_index >= 0){
									closed_elements[element_2D_index] = true;

								}

								element_2D_index = -1;

								if(prolongation_line_queue_2D.empty() && prolongation_line_queue_3D.empty()){
									break;
								}

							}
							END_TIMER("Prolongation queue");}
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
				END_TIMER("Bounding box element iteration");}

			}
		}

		END_TIMER("Element iteration");}

};


void InspectElements::computeIntersections2d3dUseProlongationTable(std::vector< unsigned int >& prolongation_table, 
                                                                   const ElementFullIter& ele_2d, 
                                                                   const ElementFullIter& ele_3d){

	for(unsigned int i = 0; i < prolongation_table.size();i++){

		unsigned int side;
		bool is_triangle_side = true;

		if(prolongation_table[i] >= 4){
			side = prolongation_table[i] - 4;
		}else{
			side = prolongation_table[i];
			is_triangle_side = false;
		}

        DBGMSG("prolongation table: %d %d\n", side, is_triangle_side);

		if(is_triangle_side){
			// prolongation through the triangle side

			//xprintf(Msg,"Procházím hranu(%d) na id elementu(%d), hrana(%d)\n"
			//					, side,elm->index(),elm->side(2-side)->el_idx());


// 			SideIter elm_side = elm->side((3-side)%3);
            SideIter elm_side = ele_2d->side(side);


			Edge *edg = elm_side->edge();

			for(int j=0; j < edg->n_sides;j++) {
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					unsigned int sousedni_element = other_side->element()->index(); // 2D element

                    DBGMSG("2d sousedni_element %d\n", sousedni_element);
							//xprintf(Msg, "Naleznut sousedni element elementu(3030) s ctyrstenem(%d)\n", ele->index());
							//xprintf(Msg, "\t\t Idx původního elementu a jeho hrany(%d,%d) - Idx sousedního elementu a jeho hrany(%d,%d)\n",elm->index(),side,other_side->element()->index(),other_side->el_idx());


						if(!intersectionExists(sousedni_element,ele_3d->index())){
							//flag_for_3D_elements[ele->index()] = sousedni_element;
                            DBGMSG("2d prolong\n");
							// Vytvoření průniku bez potřeby počítání
							IntersectionPolygon il_other(sousedni_element, ele_3d->index());
							intersection_list[sousedni_element].push_back(il_other);

							//ProlongationLine pl2(sousedni_element, elm_2D, intersection_list[sousedni_element].size() - 1);
							ProlongationLine pl2 = {sousedni_element, ele_3d->index(), intersection_list[sousedni_element].size() - 1, -1, -1};
							prolongation_line_queue_2D.push(pl2);
						}
				}
			}

		}else{
			// prolongation through the tetrahedron side

// 			SideIter elm_side = ele->side((unsigned int)(3-side)); //ele->side(3-side);
            SideIter elm_side = ele_3d->side(side);
			Edge *edg = elm_side->edge();

			for(int j=0; j < edg->n_sides;j++) {
				SideIter other_side=edg->side(j);
				if (other_side != elm_side) {
					//

					unsigned int sousedni_element = other_side->element()->index();

                    DBGMSG("3d sousedni_element %d\n", sousedni_element);

                    // TODO: flag_for_3D_elements seems to be optimalisation:
                    // - rename it, describe it and test that it is really useful !!
                    // how it probably works: 
                    // - if "-1" then no intersection has been computed for the 3D element
                    // - if not "-1" check the index of actual 2D element 
                    //   (most probable case, that we are looking along 2D element back to 3D element, which we have just computed)
                    // - if it is another 2D element, then go through all found intersections of the 3D element and test it..
                    
					if(flag_for_3D_elements[sousedni_element] == -1 || 
                        (flag_for_3D_elements[sousedni_element] != (int)ele_2d->index() && !intersectionExists(ele_2d->index(),sousedni_element))){
						
                        flag_for_3D_elements[sousedni_element] = ele_2d->index();
						// Jedná se o vnitřní čtyřstěny v trojúhelníku

                        DBGMSG("3d prolong\n");
                        
						// Vytvoření průniku bez potřeby počítání
						IntersectionPolygon il_other(ele_2d.index(), sousedni_element);
						intersection_list[ele_2d.index()].push_back(il_other);

						//ProlongationLine pl2(elm.index(), sousedni_element, intersection_list[elm.index()].size() - 1);
						ProlongationLine pl2 = {ele_2d.index(), sousedni_element, intersection_list[ele_2d.index()].size() - 1, -1, -1};
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

double InspectElements::polygon_area()
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_list.size(); i++){
        for(unsigned int j = 0; j < intersection_list[i].size();j++){
            double t2dArea = mesh->element(intersection_list[i][j].ele_2d_idx())->measure();
            double localArea = intersection_list[i][j].get_area();
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
            ElementFullIter ele = mesh->element(intersection_line_list[i][j].ele_1d_idx());            
            double t1d_length = ele->measure();
            double local_length = intersection_line_list[i][j].compute_length();
            
            arma::vec3 from = intersection_line_list[i][j][0].coords(ele);
            arma::vec3 to = intersection_line_list[i][j][1].coords(ele);
            DBGMSG("sublength from [%f %f %f] to [%f %f %f] = %f\n",
                   from[0], from[1], from[2], 
                   to[0], to[1], to[2],
                   local_length*t1d_length);
            
            subtotal += local_length*t1d_length;
        }
    }
    return subtotal;
}


} // END namespace
