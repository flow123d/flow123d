/*
 * IntersectionLocal.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include "system/system.hh"

using namespace std;
namespace computeintersection{

//class ProlongationPoint;

class IntersectionPoint{

	std::vector<double> local_coord1; // vektor lokálních souřadnic 3D elementu
	std::vector<double> local_coord2; // vektor lokální souřadnice 2D elementu

	/*	type - 2D-3D
	 *  0 = typ neurčen.
	 *  1 = SS (Stěna(3D) - Stěna(3D))
	 *  2 = SH (Stěna(3D) - Hrana(2D))
	 *  3 = HH (Hrana(2D) - Hrana(2D))
	 */
	int type;

public:
	IntersectionPoint(const std::vector<double> &c1, const double &c2, int typ = 0)
			: local_coord1(c1), local_coord2(c2), type(typ) {}

	inline std::vector<double> &el1_coord(){return local_coord1;}
	inline std::vector<double> &el2_coord(){return local_coord2;}
	//bool operator ==(const IntersectionPoint&);
};


class IntersectionLocal {

	static int numberInstance;

	int id;

	std::vector<IntersectionPoint *> i_points; //vektor ukazatelu na dvojice lokal. souradnic

	//IntersectionType type;

	unsigned int element_2D_idx;
	unsigned int element_3D_idx;

	int generateId();
public:
    typedef enum {
        point,
        line,
        area
    } IntersectionType;

    IntersectionLocal();
    IntersectionLocal(unsigned int elem2D,unsigned int elem3D);
    //Intersection_Local(IntersectionType i_type);
    //Intersection_Local(IntersectionLocal*);
    ~IntersectionLocal();

    void add_local_coord(const std::vector<double> &coordin1, const double &coordin2); //metoda na pridani souradnic do i_points
    void add_local_point(IntersectionPoint *InPoint);
       //static int getNumInstances() {
	//	return IntersectionLocal::numberInstance;
	//}

   // inline IntersectionType get_type()
    //    {return type; }
    inline IntersectionPoint * get_point(const unsigned int index)
    {
          if (index >= i_points.size() ) return NULL;
          else return i_points[index];
    }
    inline int getID(){
    	return id;
    }
    inline unsigned int idx_1D(){return element_2D_idx;}
    inline unsigned int idx_3D(){return element_3D_idx;}
    inline void print(){
    	xprintf(Msg, "ID 2D: %d ID 3D %d \n", element_2D_idx, element_3D_idx);
    	xprintf(Msg, "Local Coords: 1D: ", &i_points[0]->el2_coord());

    }
    //inline std::vector<int> getSide(const int number){
    //	if (number >= sides.size() ) return NULL;
    //	else return sides[index];
    //}

};

} // namespace computeintersection close
