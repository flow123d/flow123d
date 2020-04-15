#include "plucker.hh"

using namespace std;

Plucker::Plucker()
: coordinates_({0,0,0,0,0,0}),
  scale_(0),
  computed_(false),
  points_(3, 1, 2)
{}


Plucker::Plucker(Point a, Point b)
: Plucker()
{
    points_.set(0) = a;
    points_.set(1) = b;
    coordinates_(arma::span(0,2)) = point(1) - point(0);
    
    // Check empty
    ASSERT_DBG(arma::norm(coordinates_(arma::span(0,2)),2) > 0);

    scale_ = 0;
    scale_ = std::max(  scale_, std::fabs(coordinates_[0]));
    scale_ = std::max(  scale_, std::fabs(coordinates_[1]));
    scale_ = std::max(  scale_, std::fabs(coordinates_[2]));
    
    computed_ = false;
}

Plucker::Plucker(Point a, Point b, bool compute_pc)
: Plucker(a,b)
{
    if(compute_pc) compute();
}


double Plucker::operator*(const Plucker &b){
	return (coordinates_[0]*b[3]) + (coordinates_[1]*b[4]) + (coordinates_[2]*b[5]) + (coordinates_[3]*b[0]) + (coordinates_[4]*b[1]) +(coordinates_[5]*b[2]);
}


void Plucker::compute(){
    if(computed_) return;
    
    coordinates_[3] = coordinates_[1]*point(0)[2] - coordinates_[2]*point(0)[1];
    coordinates_[4] = coordinates_[2]*point(0)[0] - coordinates_[0]*point(0)[2];
    coordinates_[5] = coordinates_[0]*point(0)[1] - coordinates_[1]*point(0)[0];
    computed_ = true;
}

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


