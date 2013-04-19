#include "hyperplane.h"
#include <iostream>
using namespace std;




/*
HyperPlane<1,3>::HyperPlane(SPoint<3> point_a,SPoint<3> point_b):plucker(Plucker(Vector<3>(point_a,point_b),Vector<3>(point_a,point_b)*(Vector<3>)point_a)){
	  alfa = point_a;
	  beta = point_b;
	  };
*/



void HyperPlane<1,3>::toString(){
		  cout << "hyperplane<1,3>: ";
		  plucker.toString();
	  };
