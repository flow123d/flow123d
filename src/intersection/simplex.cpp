#include "simplex.h"
#include <iostream>
#include <armadillo>
using namespace std;
namespace computeintersection {

template<> Simplex<1> &Simplex<3>::getAbscissa(unsigned int index) {

	if (index < 3) {
		return Simplices[0][index];
	} else if (index < 5) {
		return Simplices[1][index - 2];
	} else {
		return Simplices[2][2];
	}

}

template<> Simplex<1> &Simplex<2>::getAbscissa(unsigned int index) {
	return Simplices[index];
}
/*template<> Simplex<1> &Simplex<1>::getAbscissa(unsigned int index) {
	return this[0];
}*/


/*template<int N> void Simplex<N,3>::setSimplex(
 Simplex<N - 1, 3> Simplex_n[N + 1]) {
 for (int i = 0; i < N + 1; i++) {
 Simplices[i] = Simplex_n[i];
 }
 }
 ;

 template <int N> void Simplex<N, 3>::toString() {
 for (int i = 0; i < N + 1; i++) {
 cout << "Simplex<" << N << ">";
 Simplices[i].toString();
 }
 }
 ;*/
/*
 Simplex<1> getAbsicca(const Simplex<2> &abs, int i){
 return abs[i];
 };

 Simplex<1> getAbsicca(const Simplex<3> &abs, int i){
 if(i < 3){
 return abs[0][i];
 }else if(i < 5){
 return abs[1][i - 2];
 }else{
 return abs[2][2];
 }
 };
 */

} // END namespace_close
