#include <iostream>
#include "spoint.h"
using namespace std;

template<int N> class Vector {};
template<> 
class Vector<3>{
private:
	double coordinates[3];
public:
	inline Vector(){};
	Vector(SPoint<3> a,SPoint<3> b);
	Vector(const double &a,const double &b,const double &c);
	explicit Vector(SPoint<3> a);
	inline ~Vector(){};
	inline double operator[](int index){return coordinates[index];};
	Vector<3> operator*(Vector<3> b);
	inline void setCoordinates(int index, double value){coordinates[index] = value;};
	void vector_product(Vector<3> a, SPoint<3> b);
	double scalar_product(Vector<3> a);
	void toString();
};
