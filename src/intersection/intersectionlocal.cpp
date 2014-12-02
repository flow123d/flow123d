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
		for(unsigned int j = 0; j < 3;j++){
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

void IntersectionLocal::traceGenericPolygon(){

	if(!is_patological){
		xprintf(Msg,"Tracing opt polygon(%d)", i_points.size());

		//xprintf(Msg,"Min: %f\n",std::numeric_limits<double>::epsilon());



		this->tracePolygonOpt();
		return;
	}else{
		xprintf(Msg,"Tracing traceConvexHull polygon(%d)", i_points.size());
		this->traceConvexHull();
		return;
	}
	xprintf(Msg,"Tracing generic polygon(%d)", i_points.size());


	std::vector<int> pp;
	std::vector<IntersectionPoint<2,3>> new_points, new_points2;

	double startA,startB,startC = 0.0;
	arma::vec::fixed<3> min; min.zeros();
	int min_index = -1;

	min[2] = 1;
	while(true){
		for(unsigned int j = 0; j < i_points.size();j++){
			if(i_points[j].getLocalCoords1()[1] >= startB && i_points[j].getLocalCoords1()[2] <= min[2]){
				if(i_points[j].getLocalCoords1()[2] == min[2]){
					if(i_points[j].getLocalCoords1()[1] <= min[1]){
						min = i_points[j].getLocalCoords1();
						min_index = j;
					}
				}else{
					min = i_points[j].getLocalCoords1();
					min_index = j;
				}
			}
		}
		if(min_index != -1){
			new_points.push_back(i_points[min_index]);
			startA = i_points[min_index].getLocalCoords1()[0];
			startB = i_points[min_index].getLocalCoords1()[1];
			startC = i_points[min_index].getLocalCoords1()[2];
			i_points.erase(i_points.begin() + min_index);
			min.zeros();min[2] = 1;
			min_index = -1;
		}else{
			break;
		}
	}

	// Po šikmé hraně
	min.zeros();min[0] = 1;
	while(true){
		for(unsigned int j = 0; j < i_points.size();j++){
			if(i_points[j].getLocalCoords1()[2] >= startC && i_points[j].getLocalCoords1()[0] <= min[0]){
				if(i_points[j].getLocalCoords1()[0] == min[0]){
					if(i_points[j].getLocalCoords1()[2] <= min[2]){
						min = i_points[j].getLocalCoords1();
						min_index = j;
					}
				}else{
					min = i_points[j].getLocalCoords1();
					min_index = j;
				}
			}
		}
		if(min_index != -1){
			new_points.push_back(i_points[min_index]);
			startA = i_points[min_index].getLocalCoords1()[0];
			startB = i_points[min_index].getLocalCoords1()[1];
			startC = i_points[min_index].getLocalCoords1()[2];
			i_points.erase(i_points.begin() + min_index);
			min.zeros();min[0] = 1;
			min_index = -1;
		}else{
			break;
		}
	}
	// dolu po CA->
	min.zeros();min[1] = 1;
	while(true){
		for(unsigned int j = 0; j < i_points.size();j++){
			if(i_points[j].getLocalCoords1()[0] >= startA && i_points[j].getLocalCoords1()[1] <= min[1]){
				if(i_points[j].getLocalCoords1()[1] == min[1]){
					if(i_points[j].getLocalCoords1()[0] <= min[0]){
						min = i_points[j].getLocalCoords1();
						min_index = j;
					}
				}else{
					min = i_points[j].getLocalCoords1();
					min_index = j;
				}
			}
		}
		if(min_index != -1){
			new_points.push_back(i_points[min_index]);
			startA = i_points[min_index].getLocalCoords1()[0];
			startB = i_points[min_index].getLocalCoords1()[1];
			startC = i_points[min_index].getLocalCoords1()[2];
			i_points.erase(i_points.begin() + min_index);
			min.zeros();min[1] = 1;
			min_index = -1;
		}else{
			break;
		}
	}

		if(new_points.size() > 1){

			new_points2.push_back(new_points[0]);

			double old1 = new_points[0].getLocalCoords1()[0];
			double old2 = new_points[0].getLocalCoords1()[1];

			for(unsigned int j = 1; j < new_points.size();j++){

				if(new_points[j].getLocalCoords1()[0] != old1 ||
						new_points[j].getLocalCoords1()[1] != old2){
					new_points2.push_back(new_points[j]);
				}
				old1 = new_points[j].getLocalCoords1()[0];
				old2 = new_points[j].getLocalCoords1()[1];
			}

			if(new_points2.size() > 1){
				if(new_points2[new_points2.size()-1].getLocalCoords1()[0] == new_points2[0].getLocalCoords1()[0] &&
					new_points2[new_points2.size()-1].getLocalCoords1()[1] == new_points2[0].getLocalCoords1()[1]){
					new_points2.pop_back();
				}
			}

		}else if(new_points.size() == 1){
			new_points2 = new_points;
		}

		i_points = new_points2;



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
			}else{
				index1 = i_points[m].getSide2();
			}

			if(i_points[n].isVertex()){
				index2 = 4 + ((i_points[n].getLocalCoords1()[0] == 1) ? 0 : ((i_points[n].getLocalCoords1()[1] == 1) ? 1 : 2));
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

void IntersectionLocal::tracePolygonOpt(){

	//return;
	fillTracingTable();

	std::vector<IntersectionPoint<2,3>> new_points;
	unsigned int start = -1;
	unsigned int end = -1;
	bool vrchol = false;

	for(unsigned int i = 0; i < 7; i++){
		if(tracing_table(i,0) != -1){
			end = i;
			start = tracing_table(i,0);
			if(tracing_table(i,1) != -1){
				new_points.push_back(i_points[tracing_table(i,1)]);
				if(i > 3){vrchol = true;}
			}
			if(!vrchol && tracing_table(i,2) != -1){
				new_points.push_back(i_points[tracing_table(i,2)]);
			}
			break;
		}
	}



	while(start != end){
		vrchol = false;
		if(tracing_table(start,1) != -1){
			new_points.push_back(i_points[tracing_table(start,1)]);
			// Vrcholy brát pouze jednou:
			if(start > 3){vrchol = true;}
		}
		if(!vrchol && tracing_table(start,2) != -1){
			new_points.push_back(i_points[tracing_table(start,2)]);
		}
		start = tracing_table(start, 0);
		//if(start == -1){return;}
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
