#include "plucker.h"

using namespace std;
namespace computeintersection{

Plucker::Plucker(){
	computed = false;
};

Plucker::Plucker(const arma::vec3 &a,const arma::vec3 &b){
	this->compute(a, b);
	computed = true;
};

Plucker::Plucker(const Plucker &p){
	coordinates = p.get_plucker_coords();
	computed = p.is_computed();
};

double Plucker::operator*(const Plucker &b){
	return (coordinates[0]*b[3]) + (coordinates[1]*b[4]) + (coordinates[2]*b[5]) + (coordinates[3]*b[0]) + (coordinates[4]*b[1]) +(coordinates[5]*b[2]);
};

void Plucker::compute(const arma::vec3 &a,const arma::vec3 &b){
	
	coordinates[0] = b[0] - a[0];
	coordinates[1] = b[1] - a[1];
	coordinates[2] = b[2] - a[2];
	coordinates[3] = coordinates[1]*a[2] - coordinates[2]*a[1];
	coordinates[4] = coordinates[2]*a[0] - coordinates[0]*a[2];
	coordinates[5] = coordinates[0]*a[1] - coordinates[1]*a[0];
	computed = true;

};

void Plucker::toString(){
	if(computed){
		cout <<"(" << coordinates[0] << "," << coordinates[1] << "," << coordinates[2] << "," << coordinates[3] << "," << coordinates[4] << "," << coordinates[5] << ")" << endl;
	}else{
		cout << "NULL ( computed = false)" << endl;
	}
};

} // END namespace_close
