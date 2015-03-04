#include <armadillo>

using namespace std;
namespace computeintersection {

/**
 * Simplex<N> represents N-dimensional simplex,
 * Simplex<0> = pointer to 3D point
 * Simplex<1> = abscissa
 * Simplex<2> = triangle
 * Simplex<3> = tetrahedron
 * Sub - simplices are made in lexicographical order
 *
 * 				Simplex<3> with 4 points 0,1,2,3 creates:
 * 								|
 * 					 -----------------------------------------------------------------------------------------------------------------------------------------------
 * 					|												|												|												|
 * 					S<2>(0,1,2) 									S<2>(0,1,3) 									S<2>(0,2,3) 									S<2>(1,2,3)
 * 						|												|												|												|
 * 		 -------------------------------				 -------------------------------				 -------------------------------				 -------------------------------
 * 		|				|				|				|				|				|				|				|				|				|				|				|
 * 		S<1>(0,1) 		S<1>(0,2) 		S<1>(1,2) 		S<1>(0,1) 		S<1>(0,3) 		S<1>(1,3) 		S<1>(0,2) 		S<1>(0,3) 		S<1>(2,3) 		S<1>(1,2) 		S<1>(1,3) 		S<1>(2,3)
 * 			|				|				|				|				|				|				|				|				|				|				|				|
 * 		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------		 -------
 * 		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|		|
 *		S<0>(0) S<0>(1) S<0>(0) S<0>(2) S<0>(1) S<0>(2) S<0>(0) S<0>(1) S<0>(0) S<0>(3) S<0>(1) S<0>(3) S<0>(0) S<0>(2) S<0>(0) S<0>(3) S<0>(2) S<0>(3) S<0>(1) S<0>(2) S<0>(1) S<0>(3) S<0>(2) S<0>(3)
 *
 * Simplex<0> has pointer to 3D point, because 3D point will be only once in memory
 *
 */
template<int N> class Simplex;

template<> class Simplex<0> {
private:
	//arma::vec3 coordinates;
	arma::vec3* coords;
public:
	inline Simplex(){
		coords = NULL;
	};
	inline Simplex(arma::vec3 **field){coords = field[0];};
	inline ~Simplex(){};

	inline void setSimplices(arma::vec3 **field){
		coords = field[0];
	};

	inline arma::vec3 &getPointCoordinates(){
		return *coords;//inates;
	};

	inline void setPointCoordinates(arma::vec3 &point){
		if(coords == NULL){
			coords = &point;
		}else{
			(*coords)[0] = point[0];
			(*coords)[1] = point[1];
			(*coords)[2] = point[2];
		}
	}

	inline void to_string(){
		cout << "Simplex<0>(" << (*coords)[0] << "," << (*coords)[1] << ","	<< (*coords)[2] << ")" << endl;
	}
};

template<int N> class Simplex {
private:
	Simplex<N - 1> Simplices[N + 1];
public:
	inline Simplex(){};

	inline Simplex(arma::vec3 **field_of_pointers_to_coordinates){

		this->setSimplices(field_of_pointers_to_coordinates);
	};

	inline void setSimplices(arma::vec3 **field_of_pointers_to_coordinates){
		arma::vec3 *temporary_pointers[N];

		for (int i = 0; i < N; i++) {
			temporary_pointers[i] = field_of_pointers_to_coordinates[i];
		};
		Simplices[0].setSimplices(temporary_pointers);
		//Simplices[0] = Simplex<N - 1>(temporary_pointers);
		for (int i = 1; i < N + 1; i++) {
			temporary_pointers[N - i] = field_of_pointers_to_coordinates[N - i + 1];
			//Simplices[i] = Simplex<N - 1>(temporary_pointers);
			Simplices[i].setSimplices(temporary_pointers);
		}
	}

	inline ~Simplex(){};

	inline Simplex<N - 1> &operator[](int a){
		return Simplices[a];
	};

	inline void to_string(){
		cout << "Simplex<" << N << ">:" << endl;
		for (int i = 0; i < N + 1; i++) {
			for (int j = 3; N <= j; j--) {
				cout << "  ";
			}
			Simplices[i].to_string();
		}
	}

	Simplex<1> &getAbscissa(unsigned int index);
};


} // END namespace_close
