#include <armadillo>

using namespace std;
namespace computeintersection{
//template<int N> class Simplex{};
//template<> class Simplex{};

template<int N> class Simplex{
private:
	Simplex<N-1> Simplices[N+1];
public:
	inline Simplex(){};
	//Simplex(arma::vec3 *field_of_coordinates);
	inline Simplex(arma::vec3 *field_of_coordinates){

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
	};
	inline ~Simplex(){};
	inline void setSimplex(Simplex<N-1> Simplex_n[N+1]){
		for (int i = 0; i < N + 1; i++) {
			Simplices[i] = Simplex_n[i];
		}
	};
	inline void setSimplex(Simplex<N-1> Simplex0){Simplices[0] = Simplex0;};
	//inline Simplex<N-1,3> operator[](int a){return Simplices[a];};
	inline Simplex<N-1> operator[](int a) const {return Simplices[a];};
	inline Simplex<N-1> getSimplexChild(int a) const {return Simplices[a];};

	inline void toString(){
		cout << "Simplex<" << N << ">:" << endl;
		for (int i = 0; i < N + 1; i++) {
			for (int j = 3; N <= j; j--){
				cout << "  ";
			}
			Simplices[i].toString();
		}};
	inline int dim(){return N;};
	inline Simplex<1> *getAbscissa(unsigned int index){

		// Nástřel
		/*
		if(N == 3){
			if(index < 3){
				return Simplices[0][index];
			}else if(index < 5){
				return &Simplices[1][index-2];
			}else{
				return &Simplices[2][2];
			}
		}else if(N == 2){
			return &Simplices[index];
		}else if(N == 1){
			return NULL;
		}else{
			return NULL;
		}*/
		return NULL;

	};
};

template<> class Simplex<0>{
private:
	arma::vec3 coordinates;
	//SPoint<3> point;
public:
	inline Simplex(){};
	inline Simplex(arma::vec3 co):coordinates(co){};
	inline Simplex(arma::vec3 *field){
		coordinates = field[0];
	};
	//inline Simplex(SPoint<3> *point_array):point(point_array[0]){};
	inline ~Simplex(){};
	inline arma::vec3 getPointCoordinates() const {return coordinates;};


	inline int dim(){return 0;};
	inline void toString(){
		cout << "Simplex<0>("<< coordinates[0] << "," << coordinates[1] << "," << coordinates[2] <<")" << endl;
	};
};
} // END namespace_close
