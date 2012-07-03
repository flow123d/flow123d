#ifndef POINT_HH_
#define POINT_HH_

#include <armadillo>

typedef arma::vec3 Vec3;

/**
 *  A point in 3d ambient space represented by an armadillo vector.
 * 
 *  TODO: implement it like  in mesh/node.*
 */ 
class Point {
public:
private:
    Vec3        point_;         /// Vector of coordinates in 3d Cartesian space.
};  

#endif