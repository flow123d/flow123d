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
#include "mesh/mesh.h"
#include <queue>


using namespace std;
namespace computeintersection{

/**
 * Class doc.
 * Naming convention.
 *
 * Rename IntersectionPolygon use only for output (traced) polygon.
 *
 * Have distinguished class for incomplete intersection polygon
 * (How it is related to ProlongationLine?)
 *
 *
 *
 */

class IntersectionLine{

	std::vector<IntersectionPoint<1,3>> i_points;
	unsigned int element_1D_idx;
	unsigned int element_3D_idx;

public:

	inline IntersectionLine(){};
	inline IntersectionLine(unsigned int ele_1D, unsigned int ele_3D):element_1D_idx(ele_1D), element_3D_idx(ele_3D){};
	inline ~IntersectionLine(){};

	inline std::vector<IntersectionPoint<1,3>> &get_points(){
		return i_points;
	}

	inline const std::vector<IntersectionPoint<1,3>> &get_points() const{
		return i_points;
	}

	inline void trace_line(){

		if(i_points.size() > 1){
			unsigned int j = (i_points[0].getSide2() + i_points[0].getOrientation())%2;

			if(j == 0){
				std::vector<IntersectionPoint<1,3>> new_points(2);
				new_points[0] = i_points[1];
				new_points[1] = i_points[0];
				i_points = new_points;
			}
		}

	};

	inline const IntersectionPoint<1,3> &operator[](unsigned int index) const{
		return i_points[index];
	};

	inline const unsigned int size() const{
		return i_points.size();
	};

	inline void add_points(const std::vector<IntersectionPoint<1,3>> &points){
		i_points = points;
	};

	inline unsigned int get_elm1D_idx() const{
		return element_1D_idx;
	};

	inline unsigned int get_elm3D_idx() const{
		return element_3D_idx;
	};

	inline const IntersectionPoint<1,3> &get_point(const unsigned int index) const
	{
		 return i_points[index];
	};

};

class IntersectionLocal {

	static const double epsilon;

	bool is_patological;

	std::vector<IntersectionPoint<2,3>> i_points; //vektor ukazatelu na dvojice lokal. souradnic
	//arma::mat::fixed<7,4> tracing_table;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;


public:

    IntersectionLocal();
    IntersectionLocal(unsigned int elem2D,unsigned int elem3D);
    ~IntersectionLocal();

    void addIP(const IntersectionPoint<2,3> &InPoint);

    inline bool isPatological(){
    	return is_patological;
    };

    inline const IntersectionPoint<2,3> &get_point(const unsigned int index) const
    {
         return i_points[index];
    }

    inline unsigned int size(){
    	return i_points.size();
    }
    inline unsigned int idx_2D(){return element_2D_idx;}
    inline unsigned int idx_3D(){return element_3D_idx;}

    /**
     * method compute local polygon area from barycentric coordinates
     * @return double computed local area
     */
    double getArea() const;

    inline void print(){
    	xprintf(Msg, "Metoda nebyla zatím implementována!\n");
    }

    /*inline void printTracingTable(){
    	tracing_table.print();
    };*/

    void traceGenericPolygon(std::vector<unsigned int> &prolongation_table);

    void trace_polygon_opt(std::vector<unsigned int> &prolongation_table);

    void trace_polygon_convex_hull(std::vector<unsigned int> &prolongation_table);

    double convex_hull_cross(const IntersectionPoint<2,3> &O,
    		const IntersectionPoint<2,3> &A,
    		const IntersectionPoint<2,3> &B) const;

    int convex_hull_prolongation_side(const IntersectionPoint<2,3> &A,
    		const IntersectionPoint<2,3> &B) const;

    int side_content_prolongation() const;
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
    /*
     * Optimalizovanější verze
     * */

    /*
     * Vrací IntersectionPoint s prohozenými dimenzemi i daty.
     * */
    template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> flipDimension(IntersectionPoint<dim, subdim> &IP){
    		//cout << "IntersectionLocal::flipDimension<" << dim << "," << subdim << "> na <" << subdim << "," <<  dim << ">" << endl;
        	IntersectionPoint<subdim, dim> IPn(IP.get_local_coords2(), IP.get_local_coords1(), IP.getSide2(), IP.getSide1(), IP.getOrientation(), IP.isVertex(), IP.isPatological());
        	return IPn;
     };

     /*
      * Interpoluje souřadnice na elementu o dimenzi nižší
      * (1 -> 2) nebo (2 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou nižší dimenzí
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-1> &IP){

        	arma::vec::fixed<d+1> interpolovane;
        	//cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-1 << "> na <" << sd << "," <<  d << ">" << endl;
        	if(d == 3){
        		interpolovane = RefSimplex<3>::interpolate<2>(IP.get_local_coords2(), IP.getSide2());
        	}else if(d == 2){
        		interpolovane = RefSimplex<2>::interpolate<1>(IP.get_local_coords2(), IP.getSide2());
        	}else{
        		cout << "zakazany stav" << endl;
        		interpolovane.zeros();
        	}
        	IntersectionPoint<sd, d> IPn(IP.get_local_coords1(),interpolovane,IP.getSide1(), IP.getSide2(),IP.getOrientation(),IP.isVertex(), IP.isPatological());
        	return IPn;
      };

     /*
      * Interpoluje souřadnice na elementu o 2 dimenze nižší
      * (1 -> 3)
      * Templetuji dimenzí, kterou chci a vkládám IP s druhou dimenzí o 2 menší
      * */
     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-2> &IP){

             	arma::vec::fixed<d+1> interpolovane;
             	//cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-2 << "> na <" << sd << "," <<  d << ">" << endl;
             	if(d == 3){
             		interpolovane = RefSimplex<3>::interpolate<1>(IP.get_local_coords2(), IP.getSide2());
             	}else{
             		cout << "zakazany stav" << endl;
             		interpolovane.zeros();
             	}
             	IntersectionPoint<sd, d> IPn(IP.get_local_coords1(),interpolovane,IP.getSide1(), IP.getSide2(), IP.getOrientation(),IP.isVertex(), IP.isPatological());
             	return IPn;
           };
};



} // namespace computeintersection close
