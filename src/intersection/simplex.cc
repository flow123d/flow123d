#include "simplex.hh"
#include <iostream>
#include <armadillo>
#include "system/system.hh"
#include "mesh/ref_element.hh"

using namespace std;
namespace computeintersection {

template<unsigned int N> void Simplex<N>::set_simplices(arma::vec3 **field_of_pointers_to_coordinates)
{
    ASSERT_DBG(field_of_pointers_to_coordinates !=nullptr);
    arma::vec3 *temporary_pointers[N];

    // filling temporary array of size N from array of size (N+1)
    for (unsigned int i = 0; i < N; i++) {
        temporary_pointers[i] = field_of_pointers_to_coordinates[i];
    };
    // Creating sub-simplices in lexicografic order
    simplices_[0].set_simplices(temporary_pointers);
    for (unsigned int i = 1; i < N + 1; i++) {
        temporary_pointers[N - i] = field_of_pointers_to_coordinates[N - i + 1];
        simplices_[i].set_simplices(temporary_pointers);
    }
}

const static std::vector<unsigned int> face_edge = {0,0,1,2,1,2};
template<> Simplex<1> &Simplex<3>::abscissa(unsigned int idx) {
    ASSERT_DBG(idx < 6);
	/* we need the first sub-simplex for getting first three abscissas
	*  the second sub-simplex for other two abscissas
	*  the third sub-simplex for the last abscissa
	*/

    unsigned int face = RefElement<3>::interact(Interaction<2,1>(idx))[0];
    return simplices_[face][face_edge[idx]];
}


template<> Simplex<1> &Simplex<2>::abscissa(unsigned int idx) {
    ASSERT_DBG(idx < 3);
	return simplices_[idx];
}

template<> Simplex<0> &Simplex<1>::node(unsigned int idx) {
    ASSERT_DBG(idx < 2);
    return simplices_[idx];
}

template<> Simplex<0> &Simplex<2>::node(unsigned int idx) {
    ASSERT_DBG(idx < 3);
    if(idx == 2) return simplices_[1][1];
    else return simplices_[0][idx];
}

template<> Simplex<0> &Simplex<3>::node(unsigned int idx) {
    ASSERT_DBG(idx < 4);
    if(idx == 3) return simplices_[1][1][1];
    else return simplices_[0].node(idx);
}

template<> ostream& operator<< <0>(ostream& os, const Simplex< 0 >& s)
{
    ASSERT_DBG(s.coords_ != nullptr);
    os << "Simplex<0>(" << (*(s.coords_))[0] << "," << (*(s.coords_))[1] << "," << (*(s.coords_))[2] << ")";
    return os;
}

template<unsigned int N> ostream& operator<<(ostream& os, const Simplex< N >& s)
{
    os << "Simplex<" << N << ">:" << endl;
        for (unsigned int i = 0; i < N + 1; i++) {
            for (unsigned int j = 3; N <= j; j--) {
                os << "  ";
            }
            os << s.simplices_[i] << endl;
        }
    return os;
}


template class Simplex<1>;
template class Simplex<2>;
template class Simplex<3>;

template ostream& operator<< <1>(ostream &os, const Simplex<1>& s); 
template ostream& operator<< <2>(ostream &os, const Simplex<2>& s); 
template ostream& operator<< <3>(ostream &os, const Simplex<3>& s); 

} // END namespace_close
