#include "hyperplane.h"


template<int N, int M> class Simplex{};

template<> class Simplex<0,3>{
private:
	SPoint<3> point;
public:
	inline Simplex(){};
	inline Simplex(SPoint<3> _point):point(_point){};
	inline Simplex(SPoint<3> *point_array):point(point_array[0]){};
	inline ~Simplex(){};
	inline SPoint<3> getPoint(){return point;};

	void toString();
};


template<int N> class Simplex<N,3>{
private:
	Simplex<N-1,3> Simplices[N+1];
	
public:
	inline Simplex(){};
	inline Simplex(SPoint<3> *point_array){
		SPoint<3> help_point[N];
		for (int i = 0; i < N; i++) {
			help_point[i] = point_array[i];
		};
		Simplices[0] = Simplex<N - 1, 3>(help_point);
		for (int i = 1; i < N + 1; i++) {
			help_point[N - i] = point_array[N - i + 1];
			Simplices[i] = Simplex<N - 1, 3>(help_point);
		}
	};
	inline ~Simplex(){};
	inline void setSimplex(Simplex<N-1,3> Simplex_n[N+1]){
		for (int i = 0; i < N + 1; i++) {
			Simplices[i] = Simplex_n[i];
		}
	};
	inline void setSimplex(Simplex<N-1,3> Simplex0){Simplices[0] = Simplex0;};
	//inline Simplex<N-1,3> operator[](int a){return Simplices[a];};
	inline Simplex<N-1,3> &operator[](int a){return Simplices[a];};
	inline void toString(){
		for (int i = 0; i < N + 1; i++) {
			cout << "Simplex<" << N << ">";
			Simplices[i].toString();
		}};
};
