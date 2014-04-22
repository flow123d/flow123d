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

	std::vector<IntersectionPoint<2,3> *> i_points; //vektor ukazatelu na dvojice lokal. souradnic

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;


public:
    IntersectionLocal();
    IntersectionLocal(unsigned int elem2D,unsigned int elem3D);
    ~IntersectionLocal();

    void add_local_coord(const std::vector<double> &coordin1, const double &coordin2); //metoda na pridani souradnic do i_points
    void add_local_point(IntersectionPoint<2,3> *InPoint);

    inline IntersectionPoint<2,3> * get_point(const unsigned int index)
    {
          if (index >= i_points.size() ) return NULL;
          else return i_points[index];
    }
    inline int getID(){
    	return id;
    }
    inline unsigned int idx_2D(){return element_2D_idx;}
    inline unsigned int idx_3D(){return element_3D_idx;}
    inline void print(){
    	//xprintf(Msg, "ID 2D: %d ID 3D %d \n", element_2D_idx, element_3D_idx);
    	//xprintf(Msg, "Local Coords: 1D: ", &i_points[0]->el2_coord());

    }
    //inline std::vector<int> getSide(const int number){
    //	if (number >= sides.size() ) return NULL;
    //	else return sides[index];
    //}

    inline static IntersectionPoint<2,3>* convertTo23(IntersectionPoint<1,2> &ip){

    	return NULL;
    }

    inline static IntersectionPoint<2,3>* convertTo23(IntersectionPoint<2,1> &ip){
       	return NULL;
    }

    inline static IntersectionPoint<2,3>* convertTo23(IntersectionPoint<1,3> &ip){
       	return NULL;
    }

};

} // namespace computeintersection close
