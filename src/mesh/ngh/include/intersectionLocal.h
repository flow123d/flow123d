#ifndef IntersectionLocalH
#define IntersectionLocalH

#include <iostream>
#include <vector>
#include <stdio.h>

using namespace std;


//class IntersectionPoint - two vectors of local coords
class IntersectionPoint {
public:
	IntersectionPoint(const std::vector<double> &c1, const std::vector<double> &c2)
		: coord1(c1), coord2(c2) {}
	IntersectionPoint(const IntersectionPoint &LC)
		: coord1(LC.coord1), coord2(LC.coord2) {}
	inline const std::vector<double> &el1_coord() const
    		{return coord1;}
	inline const std::vector<double> &el2_coord() const
			{return coord2;}
	bool operator ==(const IntersectionPoint&);
protected:
	std::vector<double> coord1;
	std::vector<double> coord2;

	friend class IntersectionLocal;
	friend IntersectionPoint *interpolate(IntersectionPoint &A1, IntersectionPoint &A2, double t);
};
IntersectionPoint *interpolate(const IntersectionPoint &A1, const  IntersectionPoint &A2, double t);


class IntersectionLocal {
public:
    typedef enum {
        point,
        line,
        area
    } IntersectionType;

    IntersectionLocal(IntersectionType i_type);
    IntersectionLocal(IntersectionLocal*);
    ~IntersectionLocal();

    void add_local_coord(const std::vector<double> &coordin1, const std::vector<double> &coordin2); //add coords to i_points
    void add_local_point(IntersectionPoint *InPoint);
    void print(FILE *out_file);

    static int getNumInstances() {
		return IntersectionLocal::numberInstance;
	}

    inline IntersectionType get_type() const
        {return type; }
    inline unsigned int n_points() const {
    	return i_points.size();
    }
    inline const IntersectionPoint * get_point(const unsigned int index) const
    {
          if (index >= i_points.size() ) return NULL;
          else return i_points[index];
    }

private:
    static int numberInstance;
    int id;

    std::vector<IntersectionPoint *> i_points; //vektor ukazatelu na dvojice lokal. souradnic
    IntersectionType type;


    int generateId();
};

#endif /* IntersectionLocalH */
