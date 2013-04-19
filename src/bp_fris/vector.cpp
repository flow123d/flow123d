#include "vector.h"
#include <iostream>
using namespace std;

// constructor, which takes two points
Vector<3>::Vector(SPoint<3> a,SPoint<3> b){
		for(int i = 0; i < 3; i++){
		coordinates[i] = b[i] - a[i];
		}
	};

// constructor, which takes three coordinates
Vector<3>::Vector(const double &a,const double &b,const double &c){
	coordinates[0] = a;
	coordinates[1] = b;
	coordinates[2] = c;
	};

// cast from a point to a vector
// explicit
Vector<3>::Vector(SPoint<3> a){
	coordinates[0] = a[0];
	coordinates[1] = a[1];
	coordinates[2] = a[2];
	};


// returns vector multiplied by another vector
Vector<3> Vector<3>::operator*(Vector<3> b){
	return Vector<3>((coordinates[1] * b[2]) - (coordinates[2] * b[1]),(coordinates[2] * b[0]) - (coordinates[0] * b[2]),(coordinates[0] * b[1]) - (coordinates[1] * b[0]));
};
	
// calculate the vector product and returns it as a new vector
void Vector<3>::vector_product(Vector<3> a, SPoint<3> b){
	coordinates[0] = (a[1] * b[2]) - (a[2] * b[1]);
	coordinates[1] = (a[2] * b[0]) - (a[0] * b[2]);
	coordinates[2] = (a[0] * b[1]) - (a[1] * b[0]);
	};

void Vector<3>::toString(){
		cout <<"(" << coordinates[0] << "," << coordinates[1] << "," << coordinates[2] << ")" << endl;
	};
