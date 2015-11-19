#include "plucker.h"

using namespace std;
namespace computeintersection{


double Plucker::operator*(const Plucker &b){
	return (coordinates_[0]*b[3]) + (coordinates_[1]*b[4]) + (coordinates_[2]*b[5]) + (coordinates_[3]*b[0]) + (coordinates_[4]*b[1]) +(coordinates_[5]*b[2]);
};

void Plucker::compute(const arma::vec3 &a,const arma::vec3 &b){
	
	coordinates_[0] = b[0] - a[0];
	coordinates_[1] = b[1] - a[1];
	coordinates_[2] = b[2] - a[2];
	coordinates_[3] = coordinates_[1]*a[2] - coordinates_[2]*a[1];
	coordinates_[4] = coordinates_[2]*a[0] - coordinates_[0]*a[2];
	coordinates_[5] = coordinates_[0]*a[1] - coordinates_[1]*a[0];
	computed_ = true;

};

ostream& operator<<(ostream& os, const Plucker& p)
{
    if(p.computed_){
        os <<"(" << p.coordinates_[0] << "," << p.coordinates_[1] << "," << p.coordinates_[2] << "," 
           << p.coordinates_[3] << "," << p.coordinates_[4] << "," << p.coordinates_[5] << ")";
    }else{
        os << "NULL (Plucker coords have not been computed)";
    }
    return os;
}

} // END namespace_close
