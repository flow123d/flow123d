/*
 * field_model_test.cpp
 *
 *  Created on: Feb 26, 2020
 *      Author: David Flanderka
 *
 *  Tests EvalPoints, Integral classes ...
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/field_model.hh"


typedef Field<3, FieldValue<3>::Scalar > ScalarField;
typedef Field<3, FieldValue<3>::VectorFixed > VectorField;


TEST(FieldModelTest, own_model) {
	ScalarField f_scal;
	VectorField f_vec;
	Fn functor;

	FieldModel<3, FieldValue<3>::VectorFixed, Fn, FieldValueCache<double>, FieldValueCache<double>> f_product(functor, f_scal.value_cache(), f_vec.value_cache());
	f_product.cache_update();
}


// ----------------------------------------------------------------------------------------------------------------
// Pass of non-fixed size matrix as argument with fixed size.
template <uint nr, uint nc>
class MyMatrix {
public:
    typedef typename arma::mat::fixed<nr, nc> MatType;

    MyMatrix(double x) {
        uint n=1;
        for (uint i=0; i<nr; ++i)
            for (uint j=0; j<nc; ++j, ++n)
                matrix_(i,j) = n*x;
    }

    arma::mat matrix() const {
        return matrix_;
    }
private:
    MatType matrix_;
};

arma::vec3 product(typename arma::Col<double>::template fixed<1> a, arma::vec3 b) {
    return b * a;
}

TEST(FieldModelTest, armamats) {
    MyMatrix<1,1> scalar(0.5);
    MyMatrix<3,1> vector(2.0);

    auto res = product(scalar.matrix(), vector.matrix());
    std::cout << "Result = " << res;
}

