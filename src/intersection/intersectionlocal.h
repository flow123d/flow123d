/*
 * IntersectionLocal.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */
//#include <algorithm>

#ifndef INTERSECTIONLOCAL_H_
#define INTERSECTIONLOCAL_H_

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
class IntersectionLocal {

	static const double epsilon;

public:

    IntersectionLocal();
    ~IntersectionLocal();


    /*
     * Vrací IntersectionPoint s prohozenými dimenzemi i daty.
     * */
    template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> flipDimension(IntersectionPoint<dim, subdim> &IP){
    		//cout << "IntersectionLocal::flipDimension<" << dim << "," << subdim << "> na <" << subdim << "," <<  dim << ">" << endl;
        	IntersectionPoint<subdim, dim> IPn(IP.get_local_coords2(), IP.get_local_coords1(), IP.get_side2(), IP.get_side1(), IP.get_orientation(), IP.is_vertex(), IP.is_patological());
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
        		interpolovane = RefSimplex<3>::interpolate<2>(IP.get_local_coords2(), IP.get_side2());
        	}else if(d == 2){
        		interpolovane = RefSimplex<2>::interpolate<1>(IP.get_local_coords2(), IP.get_side2());
        	}else{
        		cout << "zakazany stav" << endl;
        		interpolovane.zeros();
        	}
        	IntersectionPoint<sd, d> IPn(IP.get_local_coords1(),interpolovane,IP.get_side1(), IP.get_side2(),IP.get_orientation(),IP.is_vertex(), IP.is_patological());
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
             		interpolovane = RefSimplex<3>::interpolate<1>(IP.get_local_coords2(), IP.get_side2());
             	}else{
             		cout << "zakazany stav" << endl;
             		interpolovane.zeros();
             	}
             	IntersectionPoint<sd, d> IPn(IP.get_local_coords1(),interpolovane,IP.get_side1(), IP.get_side2(), IP.get_orientation(),IP.is_vertex(), IP.is_patological());
             	return IPn;
           };
};



} // namespace computeintersection close
#endif /* INTERSECTIONLOCAL_H_ */
