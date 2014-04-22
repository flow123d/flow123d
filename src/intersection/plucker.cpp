#include "plucker.h"

using namespace std;
namespace computeintersection{

Plucker::Plucker(){
	computed = false;
};

// constructor, which takes two const vectors
Plucker::Plucker(const arma::vec3 a,const arma::vec3 b){

	this->compute(a, b);
	computed = true;

};

// specific product of two pl�cker coordinates
double Plucker::operator*(Plucker b){
	return (coordinates[0]*b[3]) + (coordinates[1]*b[4]) + (coordinates[2]*b[5]) + (coordinates[3]*b[0]) + (coordinates[4]*b[1]) +(coordinates[5]*b[2]);
	};
	
// multiplied by a constant
void Plucker::operator*(double number){
		for(int i = 0;i<6;i++){
		coordinates[i] = number*coordinates[i];
		}
	};

void Plucker::setComputed(bool comp){
	computed = comp;
};

bool Plucker::isComputed(){
	return computed;
};

void Plucker::clear(){
	computed = false;
};

void Plucker::compute(const arma::vec3 a,const arma::vec3 b){
	
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
