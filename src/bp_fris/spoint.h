#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;

template<int N> class SPoint{};

template<> class SPoint<3>{
private:
	//double coordinates[3];
	arma::vec3 coordinates;
public:

	inline SPoint(){};

	//SPoint(const double& first, const double& second, const double& third);

	inline SPoint(double first, double second, double third){
		coordinates(0) = first;
		coordinates(1) = second;
		coordinates(2) = third;
	};

	inline ~SPoint(){};

	inline double operator[](int index){return coordinates(index);};

	SPoint<3> operator+(SPoint<3> other_point);

	SPoint<3> operator*(double number);

	inline arma::vec3 getCoords(){
		return coordinates;
	};

	inline std::vector<double> getCoordsVec(){
		vector<double> pole;
		pole.push_back(coordinates(0));
		pole.push_back(coordinates(1));
		pole.push_back(coordinates(2));
		return pole;
	};

	inline void setCoords(arma::vec3 souradnice){
		coordinates = souradnice;
	};

	void toString();

	void  print(std::ostream &s) const;
};
//std::ostream & operator << ( ostream &s, const Point<3> &b);
