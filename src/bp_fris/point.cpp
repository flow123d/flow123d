#include "point.h"
#include <iostream>
using namespace std;

/*SPoint<3>::SPoint(const double& first, const double& second, const double& third){
    coordinates[0] = first;
	coordinates[1] = second;
	coordinates[2] = third;
};*/

SPoint<3>::SPoint(double first, double second, double third){
    coordinates[0] = first;
	coordinates[1] = second;
	coordinates[2] = third;
};



// adds two points to one
SPoint<3> SPoint<3>::operator+(SPoint<3> point){
	return SPoint<3>(coordinates[0] + point[0], coordinates[1] + point[1], coordinates[2] + point[2]);
};

// multiplied point numeric constant
SPoint<3> SPoint<3>::operator*(double number){
	return Point<3>(number*coordinates[0],number*coordinates[1],number*coordinates[2]);
};

void SPoint<3>::toString(){
	cout <<"(" << coordinates[0] << "," << coordinates[1] << "," << coordinates[2] << ")" << endl;
};

//void Point<3>::print(std::ostream &s) const{};
//std::ostream & operator << ( std::ostream &s, const Point<3> &b) { b.print(s); }
