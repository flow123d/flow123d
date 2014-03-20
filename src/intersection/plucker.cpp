#include "plucker.h"
#include <iostream>
using namespace std;
namespace computeintersection{
// constructor, which takes two vectors
/*Plucker::Plucker(arma::vec3 u, arma::vec3 a){
	coordinates[0] = u[0];
	coordinates[1] = u[1];
	coordinates[2] = u[2];
	coordinates[3] = a[0];
	coordinates[4] = a[1];
	coordinates[5] = a[2];
		};*/

// constructor, which takes two const vectors
Plucker::Plucker(const arma::vec3 a,const arma::vec3 b){
	coordinates[0] = b[0] - a[0];
	coordinates[1] = b[1] - a[1];
	coordinates[2] = b[2] - a[2];
	coordinates[3] = coordinates[1]*a[2] - coordinates[2]*a[1];
	coordinates[4] = coordinates[2]*a[0] - coordinates[0]*a[2];
	coordinates[5] = coordinates[0]*a[1] - coordinates[1]*a[0];
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
	
void Plucker::toString(){
		cout <<"(" << coordinates[0] << "," << coordinates[1] << "," << coordinates[2] << "," << coordinates[3] << "," << coordinates[4] << "," << coordinates[5] << ")" << endl;
	};



} // END namespace_close
