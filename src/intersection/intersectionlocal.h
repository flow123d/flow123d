/*
 * IntersectionLocal.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */
//#include <algorithm>
#include "intersectionpoint.h"
#include "prolongationline.h"
#include "system/system.hh"
#include <queue>


using namespace std;
namespace computeintersection{

class IntersectionLocal {

	//static int numberInstance;
	//int id;
	bool is_patological;

	std::vector<IntersectionPoint<2,3>> i_points; //vektor ukazatelu na dvojice lokal. souradnic
	arma::mat::fixed<7,4> tracing_table;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;


public:
    IntersectionLocal();
    IntersectionLocal(unsigned int elem2D,unsigned int elem3D);
    ~IntersectionLocal();

    void addIP(IntersectionPoint<2,3> InPoint);

    inline bool isPatological(){
    	return is_patological;
    };

    inline IntersectionPoint<2,3> get_point(const unsigned int index)
    {
         return i_points[index];
    }
    //inline int getID(){
    	//return id;
    //}
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

    void traceGenericPolygon();

    void prolongationType(const IntersectionPoint<2,3> &a, const IntersectionPoint<2,3> &b, unsigned int &type, unsigned int &index) const;

    /**
     * Trasování Polygonu
     *  - po té, co se naplní trasovací tabulka se tato tabulka prochází
     *  a propojují se návaznosti.
     *
     *  Tabulka obsahuje 4 stěny a 3 vrcholy = 7 řádků.
     *  každý řádek má údaj o tom, do kterého řádku má pokračovat a
     *  zda-li je na něm 0 - 2 průniků.
     */
    void tracePolygonOpt();

    void prolongatePolygon(std::queue<ProlongationLine> &fronta2D, std::queue<ProlongationLine> &fronta3D);

    void traceConvexHull();

    double ConvexHullCross(const IntersectionPoint<2,3> &O,
    		const IntersectionPoint<2,3> &A,
    		const IntersectionPoint<2,3> &B) const;
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
    void fillTracingTable();

    /*
     * Optimalizovanější verze
     * */

    void fillTracingTable2();

    /*
     * Vrací IntersectionPoint s prohozenými dimenzemi i daty.
     * */
    template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> flipDimension(IntersectionPoint<dim, subdim> IP){
    		//cout << "IntersectionLocal::flipDimension<" << dim << "," << subdim << "> na <" << subdim << "," <<  dim << ">" << endl;
        	IntersectionPoint<subdim, dim> IPn(IP.getLocalCoords2(), IP.getLocalCoords1(), IP.getSide2(), IP.getSide1(), IP.getOrientation(), IP.isVertex(), IP.isPatological());
        	return IPn;
     };

     /*
      * Interpoluje souřadnice na elementu o dimenzi nižší
      * (1 -> 2) nebo (2 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou nižší dimenzí
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-1> IP){

        	arma::vec::fixed<d+1> interpolovane;
        	//cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-1 << "> na <" << sd << "," <<  d << ">" << endl;
        	if(d == 3){
        		interpolovane = RefSimplex<3>::interpolate<2>(IP.getLocalCoords2(), IP.getSide2());
        	}else if(d == 2){
        		interpolovane = RefSimplex<2>::interpolate<1>(IP.getLocalCoords2(), IP.getSide2());
        	}else{
        		cout << "zakazany stav" << endl;
        		interpolovane.zeros();
        	}
        	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2(),IP.getOrientation(),IP.isVertex(), IP.isPatological());
        	return IPn;
      };

     /*
      * Interpoluje souřadnice na elementu o 2 dimenze nižší
      * (1 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou dimenzí o 2 menší
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-2> IP){

             	arma::vec::fixed<d+1> interpolovane;
             	//cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-2 << "> na <" << sd << "," <<  d << ">" << endl;
             	if(d == 3){
             		interpolovane = RefSimplex<3>::interpolate<1>(IP.getLocalCoords2(), IP.getSide2());
             	}else{
             		cout << "zakazany stav" << endl;
             		interpolovane.zeros();
             	}
             	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2(), IP.getOrientation(),IP.isVertex(), IP.isPatological());
             	return IPn;
           };
};



} // namespace computeintersection close
