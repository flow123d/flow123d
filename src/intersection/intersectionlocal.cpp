/*
 * IntersectionLocal.cpp
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include <algorithm>
#include "intersectionlocal.h"

namespace computeintersection{

IntersectionLocal::IntersectionLocal(){};

IntersectionLocal::IntersectionLocal(unsigned int elem2D,unsigned int elem3D):element_2D_idx(elem2D),element_3D_idx(elem3D){

	is_patological = false;

	for(unsigned int i = 0; i < 7;i++){
		for(unsigned int j = 0; j < 4;j++){
			tracing_table(i,j) = -1;
		}
	}
};

IntersectionLocal::~IntersectionLocal(){};

void IntersectionLocal::addIP(IntersectionPoint<2,3> InPoint){
	if(InPoint.isPatological()){
		is_patological = true;
	}
	i_points.push_back(InPoint);
};

void IntersectionLocal::traceGenericPolygon(std::vector<std::pair<unsigned int, unsigned int>> &prolongation_table){

	if(!is_patological){
		xprintf(Msg,"Tracing opt polygon(%d)\n", i_points.size());
		this->tracePolygonOpt(prolongation_table);

	}else{
		xprintf(Msg,"Tracing traceConvexHull polygon(%d)", i_points.size());
		this->traceConvexHull();
	}
};

void IntersectionLocal::fillTracingTable(){

	for(unsigned int i = 0; i < i_points.size();i++){
	if(!i_points[i].isPatological()){
		if(i_points[i].getSide1() != -1){
			// jedná se o průniky 1 -> 2 resp 1 -> 3
			// Tyto průniky jsou vždy po dvojicích
			unsigned int index1, index2;
			//xprintf(Msg,"Orientace(%d), hrana(%d), stena(%d), vrchol(%d),\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
			//xprintf(Msg,"Orientace(%d), hrana(%d),stena(%d), vrchol(%d)\n", i_points[i+1].getOrientation(),i_points[i+1].getSide1(),i_points[i+1].getSide2(),i_points[i+1].isVertex());

			unsigned int j = 0;
			/* Orientace přímek podle orientace stěn
			 *
			 *  Stěny  --- Přímky
			 *  0 -> 1    0,0
			 *  0 -> 2    0,1
			 *  0 -> 3    0,0
			 *  1 -> 2    1,1
			 *  1 -> 3    1,0
			 *  2 -> 3    0,0
			 *
			 *  Zajímá nás pouze první průnik:
			 *  0 = 0
			 *  0 = 0
			 *  0 = 0
			 *  1 = 1
			 *  1 = 1
			 *  2 = 0
			 *  => stena % 2 = orientace přímky
			 *  pokud ano, primka je obracene a jen si prohodime pruniky
			 *  + pokud je hrana troj. 1, také prohodíme
			 */


			if((i_points[i].getSide2()%2) == (int)i_points[i].getOrientation()){
				j = 1;
			}

			if(i_points[i].getSide1() == 1){
				j = 1 - j;
			}

			 // pro potřebu otáčet
			unsigned int m = i + j;
			unsigned int n = i + 1 - j;

			if(i_points[m].isVertex()){
				index1 = 4 + ((i_points[m].getLocalCoords1()[0] == 1) ? 0 : ((i_points[m].getLocalCoords1()[1] == 1) ? 1 : 2));
				i_points[m].setSide1(index1-4);//uloží se pouze index vrcholu
				i_points[m].setSide2(-1);
				// Uložit si o jaky vrchol se jedna
			}else{
				index1 = i_points[m].getSide2();
			}

			if(i_points[n].isVertex()){
				index2 = 4 + ((i_points[n].getLocalCoords1()[0] == 1) ? 0 : ((i_points[n].getLocalCoords1()[1] == 1) ? 1 : 2));
				i_points[n].setSide1(index2-4);//uloží se pouze index vrcholu
				i_points[n].setSide2(-1);
			}else{
				index2 = i_points[n].getSide2();
			}

			if(tracing_table(index1,0) == -1){
				tracing_table(index1,0) = index2;
			}else{
				xprintf(Msg, "PROBLEM - na stenu(%d) s indexem další stěny(%d) se chce zapsat nova stena(%d)\n",
						index1, tracing_table(index1,0),index2);
			}


			if(m == 0 || m == 1){
				// začíná se zde polygon, přitom bod musí být brán jako koncový
				tracing_table(index1,2) = m;
			}else{
				if(tracing_table(index1,1) == -1){
					tracing_table(index1,1) = m;
				}else{
					tracing_table(index1,2) = m;
				}
			}
			if(tracing_table(index2,1) == -1){
				tracing_table(index2,1) = n;
			}else{
				tracing_table(index2,2) = n;
			}

			i++;
		}else{
			// jedná se o průniky 2 -> 1
			//xprintf(Msg,"Orientace(%d), hrana(%d), stena(%d), vrchol(%d),\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
			unsigned int stena = i_points[i].getSide2();
			unsigned int index1 = RefSimplex<3>::line_sides[stena][i_points[i].getOrientation()];
			unsigned int index2 = RefSimplex<3>::line_sides[stena][1 - i_points[i].getOrientation()];
			tracing_table(index1,0) = index2;

			if(tracing_table(index1,1) == -1){
				tracing_table(index1,1) = i;
			}else{
				tracing_table(index1,2) = i;
			}


		}
	}else{
		// Patologické případy k trasování

		// Normální pat. případ => vznikl na hraně 4stěnu


		// Vrchol ve stěně čtyřstěnů -> je potřeba ověřit index stěny čtyřstěnu


		// Vrchol celého čtyřstěnu

		//xprintf(Msg,"Patologický: Orientace(%d), hrana(%d), stena(%d), vrchol(%d)\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
			//unsigned int stena = i_points[i].getSide2();
			//unsigned int index1 = RefSimplex<3>::line_sides[stena][1];
			//unsigned int index2 = RefSimplex<3>::line_sides[stena][0];
		//xprintf(Msg, "Na stěně(%d) do stěny(%d)\n", index1, index2);
		/*tracing_table(index1,0) = index2;

		if(tracing_table(index1,1) == -1){
			tracing_table(index1,1) = i;
		}else{
			tracing_table(index1,2) = i;
		}*/


	}
	}

	//tracing_table.print();


};

void IntersectionLocal::fillTracingTable2(){

	for(unsigned int i = 0; i < i_points.size();i++){

		if(i_points[i].getSide1() != -1){
			unsigned int index1, index2;
			unsigned int j = 0;
			/* Orientace přímek podle orientace stěn
			 *
			 *  Stěny  --- Přímky
			 *  0 -> 1    0,0
			 *  0 -> 2    0,1
			 *  0 -> 3    0,0
			 *  1 -> 2    1,1
			 *  1 -> 3    1,0
			 *  2 -> 3    0,0
			 *
			 *  Zajímá nás pouze první průnik:
			 *  0 = 0
			 *  0 = 0
			 *  0 = 0
			 *  1 = 1
			 *  1 = 1
			 *  2 = 0
			 *  => stena % 2 = orientace přímky
			 *  pokud ano, primka je obracene a jen si prohodime pruniky
			 *  + pokud je hrana troj. 1, také prohodíme
			 */


			if((i_points[i].getSide2()%2) == (int)i_points[i].getOrientation()){
				j = 1;
			}

			if(i_points[i].getSide1() == 1){
				j = 1 - j;
			}

			 // pro potřebu otáčet
			unsigned int m = i + j;
			unsigned int n = i + 1 - j;

			if(i_points[m].isVertex()){
				index1 = 4 + ((i_points[m].getLocalCoords1()[0] == 1) ? 0 : ((i_points[m].getLocalCoords1()[1] == 1) ? 1 : 2));
				//i_points[m].setSide1(index1-4);//uloží se pouze index vrcholu
				i_points[m].setSide2(-1);
				// Uložit si o jaky vrchol se jedna
				// přepíšu poslední vrchol
				tracing_table(index1,1) = m;
			}else{
				index1 = i_points[m].getSide2();
			}

			if(i_points[n].isVertex()){
				index2 = 4 + ((i_points[n].getLocalCoords1()[0] == 1) ? 0 : ((i_points[n].getLocalCoords1()[1] == 1) ? 1 : 2));
				//i_points[n].setSide1(index2-4);//uloží se pouze index vrcholu
				i_points[n].setSide2(-1);
				tracing_table(index2,2) = n;
			}else{
				index2 = i_points[n].getSide2();
			}

			if(tracing_table(index1,0) == -1){
				tracing_table(index1,0) = index2;
			}else{
				xprintf(Msg, "PROBLEM - na stenu(%d) s indexem další stěny(%d) se chce zapsat nova stena(%d)\n",
						index1, tracing_table(index1,0),index2);
			}


			tracing_table(index1,2) = m; // vždy výstupní bod
			tracing_table(index2,1) = n; // vždy vstupní bod
			tracing_table(index1,3) = i_points[m].getSide1();;

			i++;
		}else{
			// jedná se o průniky 2 -> 1
			//xprintf(Msg,"Orientace(%d), hrana(%d), stena(%d), vrchol(%d),\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
			unsigned int hrana = i_points[i].getSide2();
			unsigned int index1 = RefSimplex<3>::line_sides[hrana][i_points[i].getOrientation()];
			unsigned int index2 = RefSimplex<3>::line_sides[hrana][1 - i_points[i].getOrientation()];
			tracing_table(index1,0) = index2;
			tracing_table(index1,2) = i;
			tracing_table(index2,1) = i;
		}
		}

	//tracing_table.print();


};

void IntersectionLocal::prolongationType(const IntersectionPoint<2,3> &a, const IntersectionPoint<2,3> &b, unsigned int &type, unsigned int &index) const{
	/* informace index2D, index3D
	 * typ S-S    -1      int(index 3D hrany)    => bod vznikl na hraně 4stěnu
	 * typ S-H    int(index 2D hrany)     int(index 3D steny)     => regulerni prunik na stene i hrane
	 * typ H-H    int(index vrcholu)      -1     => průnik je vrchol trojuhl.
	 *
	 * */
	int indexHrana = -1;
	int indexStena = -1;
	int indexHrana2 = -1;
	int indexStena2 = -1;

	if(a.getSide1() == -1){
		// vytahnout přes RefSimplex<3> index steny
		// záleží na orientaci -> podle ní vytáhnout 1. stenu
		indexStena = RefSimplex<3>::line_sides[a.getSide2()][(a.getOrientation()+1)%2];
	}else if(a.getSide2() == -1){
		// vytahnout přes RefSimplex<2> index hrany
		// záleží na orientaci -> podle ní vytáhnout 1. hranu
		//indexHrana = RefSimplex<2>::line_sides[a.getSide1()][(a.getOrientation()+1)%2];
		indexHrana = b.getSide1();
	}else{
		indexHrana = a.getSide1();
		indexStena = a.getSide2();
	}

	if(b.getSide1() == -1){
		// vytáhnout přes RefSimplex<3> index steny
		// zálěží na orientaci -> podle ní vytáhnout 2. stenu
		indexStena2 = RefSimplex<3>::line_sides[b.getSide2()][b.getOrientation()];
	}else if(b.getSide2() == -1){
		// vytáhnout přes RefSimplex<2> index hrany
		// záleží na orientaci -> podle ní vytáhnout 2. hranu
		//indexHrana2 = RefSimplex<2>::line_sides[b.getSide1()][b.getOrientation()];
		indexHrana2 = a.getSide1();
	}else{
		indexHrana2 = b.getSide1();
		indexStena2 = b.getSide2();
	}

	if(a.getSide2() == -1 && b.getSide2() == -1){
		indexHrana = indexHrana2 = (a.getSide1() + b.getSide1())-1;
	}

	xprintf(Msg,"Původní index hran[%d %d], sten[%d %d], orientace[%d %d]\n",a.getSide1(), b.getSide1(), a.getSide2(), b.getSide2(), a.getOrientation(), b.getOrientation());
	xprintf(Msg, "INDEX hran[%d %d], sten[%d %d]\n", indexHrana, indexHrana2,indexStena, indexStena2);
	//xprintf(Msg, "nebo hran[%d %d], sten[%d %d]\n",);

	if(indexStena != -1 && (indexStena == indexStena2)){
		xprintf(Msg,"Nové - prodlužování hranou stěny čtyřstěnu\n");
	}else if(indexHrana != -1 && (indexHrana == indexHrana2)){
		xprintf(Msg,"Nové - prodlužování hranou trojúhelníku\n");
	}else{
		xprintf(Msg,"Nové - chyba - toto se nemělo stát\n");
	}

	// Podle indexů rozhodnout, zda-li výsledná hrana polygonu
	// průsečíku je hrana 2D simplexu nebo hrana ve stěně 3D simplexu


};

void IntersectionLocal::tracePolygonOpt(std::vector<std::pair<unsigned int, unsigned int>> &prolongation_table){

	//return;
	fillTracingTable2();

	//xprintf(Msg, "\n TRASOVACÍ TABULKA \n");
	//tracing_table.print();


	std::vector<IntersectionPoint<2,3>> new_points;
	int start_idx = -1;
	int start_point_idx = -1;

	// Nalezení prvního neprázdného řádku
	for(unsigned int i = 0; i < 7;i++){

		if(tracing_table(i,0) != -1){
			start_idx = i;
			start_point_idx = tracing_table(i,1);
			break;
		}

	}

	// trasování - max 7 iterací, kvůli eliminaci smyčky
	for(unsigned int i = 0; i < 7 && start_idx != -1;i++){

		if(start_idx < 4){
			// Jedná se o průsečíky na stěnách čtyřstěnu

			if(tracing_table(start_idx,1) == start_point_idx && i > 0){
				break;
			}else{
				xprintf(Msg,"\t\t přidán bod(%d) -- ", (int)tracing_table(start_idx,1));
				new_points.push_back(i_points[tracing_table(start_idx,1)]);
			}

			xprintf(Msg, "\t Prodlužuji stěnou(%d) - body[%d,%d]\n", start_idx,(int)tracing_table(start_idx,1),(int)tracing_table(start_idx,2));

			prolongation_table.push_back(std::make_pair(start_idx, 1));//IntersectionLocal::PROLONGATION_TYPE_TETRAHEDRON_SIDE));
			// Prodloužení -> podle indexu steny vratim sousedni 3D element

			/*SideIter elm_side = element_3D->side(start_idx);
			Edge *edg = elm_side->edge();

			if(edg == NULL){
				// Není už žádný soused
			}else{
				for(int j=0; j < edg->n_sides;j++) {
					SideIter other_side=edg->side(j);
					if (other_side != elm_side) {
						xprintf(Msg, "\t\t Idx původního elementu a jeho stěny(%d,%d) - Idx sousedního elementu a jeho stěny(%d,%d)",element_3D->index(),start_idx,other_side->element()->index(),other_side->el_idx());
					}
				}
			}*/



		}

		if(tracing_table(start_idx,3) != -1){
			// body jsou spojeny hranou trojúhelníka


			if(tracing_table(start_idx,2) == start_point_idx && i > 0){

				break;
			}else{
				xprintf(Msg,"\t\t přidán bod(%d) -- ", (int)tracing_table(start_idx,2));
				new_points.push_back(i_points[tracing_table(start_idx,2)]);
			}
			xprintf(Msg, "\t Prodlužuji hranou(%d) - body[%d,%d]\n", (int)tracing_table(start_idx,3),(int)tracing_table(start_idx,2),(int)tracing_table(tracing_table(start_idx,0),1));
			// Prodloužení -> podle indexu hrany vratim sousedni 2D element
			prolongation_table.push_back(std::make_pair(tracing_table(start_idx,3), 0));

			/*SideIter elm_side = element_2D->side(tracing_table(start_idx,3));
			Edge *edg = elm_side->edge();

			if(edg == NULL){
				// žádný soused
			}else{
				for(int j=0; j < edg->n_sides;j++) {
					SideIter other_side=edg->side(j);
					if (other_side != elm_side) {
						xprintf(Msg, "\t\t Idx původního elementu a jeho hrany(%d,%d) - Idx sousedního elementu a jeho hrany(%d,%d)",element_2D->index(),(int)tracing_table(start_idx,3),other_side->element()->index(),other_side->el_idx());
					}
				}
			}*/
		}


		start_idx = tracing_table(start_idx,0);

	}



	i_points = new_points;

};

void IntersectionLocal::traceConvexHull(){


	if(i_points.size() <= 1){
			return;
	}

	std::sort(i_points.begin(),i_points.end());

	// Odstranění duplicit
	for(unsigned int i = 0; i < i_points.size()-1; i++){
		if((i_points[i].getLC1()[0] == i_points[i+1].getLC1()[0]) &&
				(i_points[i].getLC1()[1] == i_points[i+1].getLC1()[1])){
			i_points.erase(i_points.begin()+i);
		}
	}

	if(i_points.size() > 1 && i_points[0].getLC1()[0] == i_points[i_points.size()-1].getLC1()[0] &&
			i_points[0].getLC1()[1] == i_points[i_points.size()-1].getLC1()[1]){
		i_points.erase(i_points.end());
	}

	int n = i_points.size(), k = 0;
	std::vector<IntersectionPoint<2,3>> H(2*n);

	for(int i = 0; i < n; ++i){
		while(k >= 2 && ConvexHullCross(H[k-2], H[k-1], i_points[i]) <= 0) k--;
		H[k++] = i_points[i];
	}

	for(int i = n-2, t = k+1;i>=0;i--){
		while(k >= t && (ConvexHullCross(H[k-2], H[k-1], i_points[i])) <= 0) k--;
		H[k++] = i_points[i];
	}

	H.resize(k);
	i_points = H;
};

double IntersectionLocal::ConvexHullCross(const IntersectionPoint<2,3> &O,
		const IntersectionPoint<2,3> &A,
		const IntersectionPoint<2,3> &B) const{
	return ((A.getLC1()[1]-O.getLC1()[1])*(B.getLC1()[2]-O.getLC1()[2])
			-(A.getLC1()[2]-O.getLC1()[2])*(B.getLC1()[1]-O.getLC1()[1]));
}

/* split the polygon into triangles according to the first point
 * for every triangle compute area from barycentric coordinates
 * 				Barycentric coords. 		Local coords.
 * Point A		[A0,A1,A2]					[A1,A2,0]
 * Point B	    [B0,B1,B2]					[B1,B2,0]
 * Point C		[C0,C1,C2]					[C1,C2,0]
 *
 *  triangle area: (B-A)x(C-A)/2 => ((B1-A1)(C2-A2) - (B2-A2)(C1-A1))/2
 *  => (A1(B2-C2) + B1(C2-A2) + C1(A2-B2))/2
 *
 *  final polygon area is sum of triangle areas
 */
double IntersectionLocal::getArea(){
	double subtotal = 0.0;
	for(unsigned int j = 2; j < i_points.size();j++){
		//xprintf(Msg, "volani %d %d\n",j, i_points.size());
		subtotal += fabs(i_points[0].getLocalCoords1()(1)*(i_points[j-1].getLocalCoords1()(2) - i_points[j].getLocalCoords1()(2)) +
				 i_points[j-1].getLocalCoords1()(1)*(i_points[j].getLocalCoords1()(2) - i_points[0].getLocalCoords1()(2)) +
				 i_points[j].getLocalCoords1()(1)*(i_points[0].getLocalCoords1()(2) - i_points[j-1].getLocalCoords1()(2)));
	}
	return fabs(subtotal/2);
};


} // namespace computeintersection close
