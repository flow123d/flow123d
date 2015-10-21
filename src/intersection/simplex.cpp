#include "simplex.h"
#include <iostream>
#include <armadillo>
using namespace std;
namespace computeintersection {

template<> Simplex<1> &Simplex<3>::getAbscissa(unsigned int index) {

	/* we need the first sub-simplex for getting first three abscissas
	*  the second sub-simplex for other two abscissas
	*  the third sub-simplex for the last abscissa
	*/
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


} // END namespace_close
