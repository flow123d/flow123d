#include <iostream>
using namespace std;

template<int N> class SPoint{};

template<> class SPoint<3>{
private:
	double coordinates[3];
public:

	inline SPoint(){};

	//SPoint(const double& first, const double& second, const double& third);

	SPoint(double first, double second, double third);

	inline ~SPoint(){};

	inline double operator[](int index){return coordinates[index];};

	SPoint<3> operator+(SPoint<3> other_point);

	SPoint<3> operator*(double number);

	void toString();

	void  print(std::ostream &s) const;
};
//std::ostream & operator << ( ostream &s, const Point<3> &b);
