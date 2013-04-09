/*
 * IntersectionLocal.h
 *
 *  Created on: 27.3.2013
 *      Author: viktor
 */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <queue>

//#include "system/system.hh"
//#include "system/file_path.hh"
#include "mesh/mesh.h"
//#include "mesh/msh_gmshreader.h"
//#include "mesh/bih_tree.hh"
//#include "fields/field_interpolated_p0.hh"
//#include "system/sys_profiler.hh"




using namespace std;

namespace fast_1_3{

class ProlongationPoint;

class IntersectionPoint{

	std::vector<double> local_coord1; // vektor lokálních souřadnic 3D elementu
	double local_coord2; // vektor lokální souřadnice 1D elementu

public:
	IntersectionPoint(const std::vector<double> &c1, const double &c2)
			: local_coord1(c1), local_coord2(c2) {}
	IntersectionPoint(const IntersectionPoint &LC)
			: local_coord1(LC.local_coord1), local_coord2(LC.local_coord2) {}

	inline std::vector<double> &el1_coord(){return local_coord1;}
	inline double &el2_coord(){return local_coord2;}
	//bool operator ==(const IntersectionPoint&);
};


class IntersectionLocal {

		static int numberInstance;

	    int id;

	    std::vector<IntersectionPoint *> i_points; //vektor ukazatelu na dvojice lokal. souradnic

	    //IntersectionType type;

	    unsigned int element_1D_idx;
	    unsigned int element_3D_idx;

	    int generateId();
public:
    typedef enum {
        point,
        line,
        area
    } IntersectionType;

    IntersectionLocal();
    IntersectionLocal(unsigned int elem1D,unsigned int elem3D);
    //Intersection_Local(IntersectionType i_type);
    //Intersection_Local(IntersectionLocal*);
    ~IntersectionLocal();

    //void set_elements(TElement *elem1, TElement *elem2); //metoda na naplneni ele1, ele2
    void add_local_coord(const std::vector<double> &coordin1, const double &coordin2); //metoda na pridani souradnic do i_points
    void add_local_point(IntersectionPoint *InPoint);
    void print(FILE *out_file);

    //void AddNewLocalcoord(); //doplnit predavany parametr, pridat novou i_points

    static int getNumInstances() {
		return IntersectionLocal::numberInstance;
	}

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
    /*inline std::vector<int> getSide(const int number){
    	if (number >= sides.size() ) return NULL;
    	else return sides[index];
    }*/

};

} // namespace fast_1_3 close


