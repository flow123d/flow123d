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
#include <armadillo>
#include <utility>
#include "arma_expect.hh"

#include "fields/field.hh"
#include "fields/field_common.hh"
#include "fields/field_values.hh"
#include "fields/field_value_cache.hh"


template<int spacedim, class Value>
class FieldCached /*: public FieldCommon*/ {
public:
	/// Constructor
	FieldCached()
	: fvc(Value::NRows_, Value::NCols_) {}

protected:
    FieldValueCache<typename Value::element_type> fvc;
};


template<int spacedim, class Value, class Fn, class ... Args>
class FieldModel : FieldCached<spacedim, Value> {
private:
    //typedef decltype(Fn(std::declval<Args...>())) FnReturnType;
    std::tuple<Args...> inputs;

public:
    FieldModel(Args... args)
    : inputs( std::make_tuple(std::forward<Args>(args)...) )
    { static_assert( std::is_same<typename Value::return_type, typename Fn::Result>::value, "Non-convertible functor type!"); }

    void cache_update(ElementCacheMap &cache_map) {
        // use tuple to functor
    }

};


using Vector = arma::vec3;
using Tensor = arma::mat33;
typedef Field<3, FieldValue<3>::Scalar > ScalarField;
typedef Field<3, FieldValue<3>::VectorFixed > VectorField;


class Fn {
public:
    typedef Vector Result;
    typedef double Param0;
    typedef Vector Param1;
    typedef std::tuple< std::shared_ptr<FieldCached<3, Param0>>, std::shared_ptr<FieldCached<3, Param1>> > DepFields;

    Result operator() (Param0 a, Param1 v) {
        return a * v;
    }
};


TEST(FieldModelTest, own_model) {
	ScalarField f_scalar;
	VectorField f_vec;

	FieldModel<3, FieldValue<3>::VectorFixed, Fn, ScalarField, VectorField> f_product(f_scalar, f_vec);
}
