/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "inspectelements.h"
namespace computeintersection {

InspectElements::InspectElements(){};

InspectElements::InspectElements(Mesh* _mesh):mesh(_mesh){
	//projeti.assign(sit->n_elements(),false);

	unsigned int elementLimit = 20;
		BIHTree bt(mesh, elementLimit);

		Profiler::initialize(MPI_COMM_WORLD);
		{ START_TIMER("Test_bez_BIHtree");

		FOR_ELEMENTS(mesh, elm) {
			// Pouze pro dvojice 2d - 3d
			if(elm->dim() == 3){
				xprintf(Msg, "-----Nalezen 3D element------ \n");

				this->UpdateTetrahedron(elm);
			}

			if (elm->dim() == 2) {
				xprintf(Msg, "-----Nalezen 2D element------ \n");
				this->UpdateTriangle(elm);
			 }
		}


		ComputeIntersection<Simplex<2>,Simplex<3> > CI_23(triangle, tetrahedron);
		CI_23.init();
		CI_23.compute(temporary_intersection);
		all_intersections.push_back(temporary_intersection);

		END_TIMER("Test_bez_BIHtree");}
		Profiler::instance()->output(cout);
		Profiler::uninitialize();

};

InspectElements::~InspectElements(){};


void InspectElements::UpdateTriangle(const ElementFullIter &element_2D){

			arma::vec3 pole_bodu[3] = {element_2D->node[0]->point(),element_2D->node[1]->point(),element_2D->node[2]->point()};
		    triangle = Simplex<2>(pole_bodu);

};

void InspectElements::UpdateTetrahedron(const ElementFullIter &element_3D){

			arma::vec3 pole_bodu[4] = {element_3D->node[0]->point(),element_3D->node[1]->point(),element_3D->node[2]->point(),element_3D->node[3]->point()};
		    tetrahedron = Simplex<3>(pole_bodu);

};

void InspectElements::print(char *nazev, unsigned int vyber = 0){

	FILE * soubor;
	soubor = fopen(nazev,"w");


	fprintf(soubor, "$MeshFormat\n");
	fprintf(soubor, "2.2 0 8\n");
	fprintf(soubor, "$EndMeshFormat\n");
	fprintf(soubor, "$Nodes\n");

	int pocet_pruseciku = 0;
		for(unsigned int i = 0; i < all_intersections.size(); i++){
			pocet_pruseciku += all_intersections[i].getIPsize();
		}
	fprintf(soubor, "%d\n", pocet_pruseciku);

	int velikost = 0;
	int pocet_prus = 0;
	for(unsigned int i = 0; i < all_intersections.size(); i++){


		IntersectionLocal il = all_intersections[i];
		velikost = il.getIPsize();

			ElementFullIter el2D = mesh->element(all_intersections[i].idx_2D());
			ElementFullIter el3D = mesh->element(all_intersections[i].idx_3D());

			for(unsigned int j = 0; j < velikost; j++){
				pocet_prus++;

				IntersectionPoint<2,3> IP23 = il.get_point(j);

				arma::vec3 T_globalni = (IP23.getLocalCoords1())[0] * el2D->node[0]->point()
													   +(IP23.getLocalCoords1())[1] * el2D->node[1]->point()
													   +(IP23.getLocalCoords1())[2] * el2D->node[2]->point();
				arma::vec3 C_globalni = (IP23.getLocalCoords2())[0] * el3D->node[0]->point()
													   +(IP23.getLocalCoords2())[1] * el3D->node[1]->point()
													   +(IP23.getLocalCoords2())[2] * el3D->node[2]->point()
													   +(IP23.getLocalCoords2())[3] * el3D->node[3]->point();

				fprintf(soubor,"%d %f %f %f\n", pocet_prus, T_globalni[0], T_globalni[1], T_globalni[2]);
				//fprintf(soubor,"%d %f %f %f\n", pocet_prus, C_globalni[0], C_globalni[1], C_globalni[2]);

			}
	}

	fprintf(soubor,"$EndNodes\n");
	fprintf(soubor,"$Elements\n");
	fprintf(soubor,"%d\n", pocet_pruseciku);

	pocet_prus = 0;
	int posledni = 0;
		for(unsigned int i = 0; i < all_intersections.size(); i++){


			IntersectionLocal il = all_intersections[i];
			velikost = il.getIPsize();

				for(unsigned int j = 0; j < velikost; j++){
					pocet_prus++;
					if(j == 0){
						posledni = pocet_prus;
					}

					if((j+1) == velikost){
						fprintf(soubor,"%d 1 2 18 7 %d %d\n", pocet_prus, pocet_prus, posledni);
					}else{
						fprintf(soubor,"%d 1 2 18 7 %d %d\n", pocet_prus, pocet_prus, pocet_prus+1);
					}

				}
		}

	fprintf(soubor,"$EndElements\n");
	fclose(soubor);

}

} // END namespace
