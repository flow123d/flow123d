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
    inline void print(){
    	//xprintf(Msg, "ID 2D: %d ID 3D %d \n", element_2D_idx, element_3D_idx);
    	//xprintf(Msg, "Local Coords: 1D: ", &i_points[0]->el2_coord());

    }



    template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> flipDimension(IntersectionPoint<dim, subdim> IP){
    		cout << "IntersectionLocal::flipDimension<" << dim << "," << subdim << "> na <" << subdim << "," <<  dim << ">" << endl;
        	IntersectionPoint<subdim, dim> IPn(IP.getLocalCoords2(), IP.getLocalCoords1(), IP.getSide2(), IP.getSide1());
        	return IPn;
     };

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
        	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2());
        	return IPn;
      };

     template<int sd,int d> inline static IntersectionPoint<sd, d> interpolateDimension(IntersectionPoint<sd,d-2> IP){

             	arma::vec::fixed<d+1> interpolovane;
             	cout << "IntersectionLocal::interpolateDimension<" << sd << "," << d-2 << "> na <" << sd << "," <<  d << ">" << endl;
             	if(d == 3){
             		interpolovane = RefSimplex<3>::interpolate<1>(IP.getLocalCoords2(), IP.getSide2());
             	}else{
             		cout << "zakazany stav" << endl;
             		interpolovane.zeros();
             	}
             	IntersectionPoint<sd, d> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2());
             	return IPn;
           };

        /*template<unsigned int subdim, unsigned int dim> inline static IntersectionPoint<subdim, dim> interpolateDimension(IntersectionPoint<subdim, dim-2> &IP){

    		arma::vec::fixed<dim+1> interpolovane = RefSimplex<dim>::interpolate<dim-2>(IP.getSide2());
    		IntersectionPoint<subdim, dim> IPn(IP.getLocalCoords1(),interpolovane,IP.getSide1(), IP.getSide2());

    		return IPn;
    	};*/
    //inline std::vector<int> getSide(const int number){
    //	if (number >= sides.size() ) return NULL;
    //	else return sides[index];
    //}





};



} // namespace computeintersection close
