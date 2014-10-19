/*
 * IntersectionLocal.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */
#include "intersectionpoint.h"
#include "system/system.hh"

using namespace std;
namespace computeintersection{

class IntersectionLocal {

	//static int numberInstance;
	int id;

	std::vector<IntersectionPoint<2,3>> i_points; //vektor ukazatelu na dvojice lokal. souradnic
	arma::mat::fixed<7,3> tracing_table;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;


public:
    IntersectionLocal();
    IntersectionLocal(unsigned int elem2D,unsigned int elem3D);
    ~IntersectionLocal();

    void addIP(IntersectionPoint<2,3> InPoint);

    inline IntersectionPoint<2,3> get_point(const unsigned int index)
    {
         return i_points[index];
    }
    inline int getID(){
    	return id;
    }
    inline int getIPsize(){
    	return i_points.size();
    }
    inline unsigned int idx_2D(){return element_2D_idx;}
    inline unsigned int idx_3D(){return element_3D_idx;}

    /**
     * method compute local polygon area from barycentric coordinates
     * @return double computed local area
     */
    double getArea();

    inline void print(){
    	xprintf(Msg, "Metoda nebyla zatím implementována!\n");
    }

    inline void printTracingTable(){
    	tracing_table.print();
    };

    inline void setTracingTable(unsigned int rows, unsigned int cols, unsigned int value){
    	tracing_table(rows, cols) = value;
    };

    inline unsigned int getTracingTableValue(unsigned int rows, unsigned int cols){
    	return tracing_table(rows, cols);
    };

    /**
     * Trasování Polygonu
     *  - po té, co se naplní trasovací tabulka se tato tabulka prochází
     *  a propojují se návaznosti.
     *
     *  Tabulka obsahuje 4 stěny a 3 vrcholy = 7 řádků.
     *  každý řádek má údaj o tom, do kterého řádku má pokračovat a
     *  zda-li je na něm 0 - 2 průniků.
     */
    inline void tracePolygon(){

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

    /*
     * Naplnění trasovací tabulky
     *  A) nejdříve se procházejí průniky od přímek trojúhelníku
     *   - zjišťuje se, jak je orientovaná přímka vůči stěnám.
     *   pokud je obráceně, body se prohodí + pokud se jedná o přímku
     *   s indexem 1. také se prohodí.
     *   - zjištuje se, jestli se jedná o průnik na stěně nebo na vrcholu
     *   - pokud se jedná o první průnik = jedná se i o koncový průnik a
     *   musí být zapsán do 2.fáze průniků.
     *
     *   - Každý průnik se zapisuje ke stěně/vrcholu, kde vznik (krom prvního),
     *   a zapisuje se k prvnímu průniku stěna/vrchol na který pokračuje
     *
     * 	B) prochází se průniky od hran čtyřstěnu
     * 	 - podle orientací průniku + orientace stěn se vybere stěna
     * 	 z referenčního elementu, ke které průnik patří a do které
     * 	 průnik pokračuje
     * */
    inline void fillTracingTable(){

    	for(unsigned int i = 0; i < i_points.size();i++){

    		if(i_points[i].getSide1() != -1){
    			// jedná se o průniky 1 -> 2 resp 1 -> 3
    			// Tyto průniky jsou vždy po dvojicích
    			unsigned int index1, index2;
    			xprintf(Msg,"Orientace(%d), hrana(%d), stena(%d), vrchol(%d),\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
    			xprintf(Msg,"Orientace(%d), hrana(%d),stena(%d), vrchol(%d)\n", i_points[i+1].getOrientation(),i_points[i+1].getSide1(),i_points[i+1].getSide2(),i_points[i+1].isVertex());

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
    			xprintf(Msg,"Orientace(%d), hrana(%d), stena(%d), vrchol(%d),\n", i_points[i].getOrientation(),i_points[i].getSide1(),i_points[i].getSide2(),i_points[i].isVertex());
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
    	}

    	tracing_table.print();
    }

    /*
     * Vrací IntersectionPoint s prohozenými dimenzemi i daty.
     * */
    template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> flipDimension(IntersectionPoint<dim, subdim> IP){
    		cout << "IntersectionLocal::flipDimension<" << dim << "," << subdim << "> na <" << subdim << "," <<  dim << ">" << endl;
        	IntersectionPoint<subdim, dim> IPn(IP.getLocalCoords2(), IP.getLocalCoords1(), IP.getSide2(), IP.getSide1(), IP.getOrientation(), IP.isVertex());
        	return IPn;
     };

     /*
      * Interpoluje souřadnice na elementu o dimenzi nižší
      * (1 -> 2) nebo (2 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou nižší dimenzí
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-1> IP){

        	arma::vec::fixed<d+1> interpolovane;
        	cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-1 << "> na <" << sd << "," <<  d << ">" << endl;
        	if(d == 3){
        		interpolovane = RefSimplex<3>::interpolate<2>(IP.getLocalCoords2(), IP.getSide2());
        	}else if(d == 2){
        		interpolovane = RefSimplex<2>::interpolate<1>(IP.getLocalCoords2(), IP.getSide2());
        	}else{
        		cout << "zakazany stav" << endl;
        		interpolovane.zeros();
        	}
        	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2(),IP.getOrientation(),IP.isVertex());
        	return IPn;
      };

     /*
      * Interpoluje souřadnice na elementu o 2 dimenze nižší
      * (1 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou dimenzí o 2 menší
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-2> IP){

             	arma::vec::fixed<d+1> interpolovane;
             	cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-2 << "> na <" << sd << "," <<  d << ">" << endl;
             	if(d == 3){
             		interpolovane = RefSimplex<3>::interpolate<1>(IP.getLocalCoords2(), IP.getSide2());
             	}else{
             		cout << "zakazany stav" << endl;
             		interpolovane.zeros();
             	}
             	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2(), IP.getOrientation(),IP.isVertex());
             	return IPn;
           };
};



} // namespace computeintersection close
