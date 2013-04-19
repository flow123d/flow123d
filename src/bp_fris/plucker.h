#include <iostream>
#include "vector.h"
using namespace std;

class Plucker{
private:
	double coordinates[6];
public:
	inline Plucker(){};
	Plucker(Vector<3> u, Vector<3> a);
	inline ~Plucker(){};
	inline double operator[](int index){return coordinates[index];};
	double operator*(Plucker b);
	void operator*(double number);
	void toString();
};
