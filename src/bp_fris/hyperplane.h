#include <iostream>
#include "plucker.h"
using namespace std;


template<int N, int M> class HyperPlane{};

template<> class HyperPlane<0,3>{
private:
	SPoint<3> point;
public:
	inline HyperPlane(){};
	inline HyperPlane(SPoint<3> _point):point(_point){};
	inline ~HyperPlane(){};
	inline SPoint<3> getPoint(){return point;};
};

template<> class HyperPlane<1,3>{
  private:
	  Plucker plucker;
	  SPoint<3> alfa;
	  SPoint<3> beta;
  public:
	  inline HyperPlane(){};
	  inline HyperPlane(SPoint<3> point_a,SPoint<3> point_b):plucker(Plucker(Vector<3>(point_a,point_b),Vector<3>(point_a,point_b)*(Vector<3>)point_a)){
		  alfa = point_a;
		  beta = point_b;
	  };
	  inline ~HyperPlane(){};
	  inline Plucker getPlucker(){return plucker;};
	  inline SPoint<3> getPointA(){return alfa;};
	  inline SPoint<3> getPointB(){return beta;};


	  void toString();
};
