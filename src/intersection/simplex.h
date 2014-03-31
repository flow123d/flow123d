#include <armadillo>

using namespace std;
namespace computeintersection {

// Declaration of template Simplex - for compiler
template<int N> class Simplex;

template<> class Simplex<0> {
private:
	arma::vec3 coordinates;
	//SPoint<3> point;
public:
	inline Simplex() {
	}
	;
	inline Simplex(arma::vec3 co) :
			coordinates(co) {
	}
	;
	inline Simplex(arma::vec3 *field) {
		coordinates = field[0];
	}
	;
	//inline Simplex(SPoint<3> *point_array):point(point_array[0]){};
	inline ~Simplex() {
	}
	;
	inline arma::vec3 getPointCoordinates() const {
		return coordinates;
	}
	;

	inline int dim() {
		return 0;
	}
	;
	inline void toString() {
		cout << "Simplex<0>(" << coordinates[0] << "," << coordinates[1] << ","
				<< coordinates[2] << ")" << endl;
	}
	;
};

template<int N> class Simplex {
private:
	Simplex<N - 1> Simplices[N + 1];
public:
	inline Simplex() {}

	//Simplex(arma::vec3 *field_of_coordinates);
	inline Simplex(arma::vec3 *field_of_coordinates) {

		arma::vec3 help_coordinates[N];
		//SPoint<3> help_point[N];
		for (int i = 0; i < N; i++) {
			help_coordinates[i] = field_of_coordinates[i];
		};
		Simplices[0] = Simplex<N - 1>(help_coordinates);
		for (int i = 1; i < N + 1; i++) {
			help_coordinates[N - i] = field_of_coordinates[N - i + 1];
			Simplices[i] = Simplex<N - 1>(help_coordinates);
		}
	}
	;
	inline ~Simplex() {
	}
	;
	inline void setSimplex(Simplex<N - 1> Simplex_n[N + 1]) {
		for (int i = 0; i < N + 1; i++) {
			Simplices[i] = Simplex_n[i];
		}
	}
	;
	inline void setSimplex(Simplex<N - 1> Simplex0) {
		Simplices[0] = Simplex0;
	}
	;
	//inline Simplex<N-1,3> operator[](int a){return Simplices[a];};
	inline Simplex<N - 1> operator[](int a) const {
		return Simplices[a];
	}
	;
	inline Simplex<N - 1> getSimplexChild(int a){
		return Simplices[a];
	}
	;

	inline void toString() {
		cout << "Simplex<" << N << ">:" << endl;
		for (int i = 0; i < N + 1; i++) {
			for (int j = 3; N <= j; j--) {
				cout << "  ";
			}
			Simplices[i].toString();
		}
	}
	;
	inline int dim() {
		return N;
	}
	;
	Simplex<1> getAbscissa(unsigned int index);
};


} // END namespace_close
